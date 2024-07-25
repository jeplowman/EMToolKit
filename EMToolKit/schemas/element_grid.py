# This element_grid class builds on coord_grid to specify a grid of basis or detector
# elements (e.g., the spatial response functions of an imaging detector, including
# their pixel boxes and PSF). It is the exemplar class for the source and detector
# objects used by get_sparse_response_matrix.
# The primary methods of these objects are as follows:
#    elements: Returns the properties of the element(s) at address i.
#              When called as [elms,vals,pnts] = some_element_grid.elements(i),
#              must return the following:
#                  pnts: The coordinate points at which the these elements are non-zero,
#                        in the element_grid's coordinate frame. Dimensions should be
#                        npts by ndim, where ndim is the dimensionality of the
#                        source coordinate system.
#                  vals: The values of the elements' 'basis' function(s) at
#                        those coordinates. Dimensions are npts.
#                  elms: Indices of the element(s) corresponding to each of those
#                        point/value pairs. Most often these will all be the same as i,
#                        but they don't have to be. Dimensions are npts.
#    response: Returns the response of the elements to a delta function source at a given
#              point in the element_grid's coordinate frame. When called as
#              [vals,elms] = detector.response(point), must return the following:
#                  elms: The indices of every element in the grid which has a non-zero
#                        response to a delta function source at the given point.
#                  vals: The values of each of those responses. Same dimensions as elms.
#
# The objects also have the following attributes which are intended for external use:
#    coords: Information about the coordinate system of the output of elements and
#            the input to response. Implemented as an instance of coordgrid.
#    nelm:   The number of elements in the element_grid, as well as the number
#            of unique element id/indices.
#    nadr:   The number of element addresses in the element grid. These are
#            how the elements are accessed via the elements method. In
#            basic use these will be the same as the element indices,
#            but there are use cases where other addressing modes are
#            desired. The addresses are generally expected to be positive integers.
# Elements are addressed by a single index apiece, but
# this is trivial to multiplex under the hood using, for instance,
# i = i0*nj*nk+j0*nk+k0. Numpy's ravel_index, unravel_multi_index, and
# flatten functions provide exactly this functionality. The intention of this class
# is to implement the functionalite needed by get_sparse_response_matrix for typical
# detector arrays e.g., 2D imagers or 3D spectrographs, as well as for gridded arrays
# of source elements. It should be fairly powerful and extensible.

# There is an issue with rounding and floating point jitter for aligned grids that
# share some ratios in their spacing. We add this small offset when discretizing
# grid indices to avoid this, which is a bit of a hack...
tiny = 1.0e-4
# On a related note, there are some intrinsic issues with even nsubgrid, which
# can produce subpixel shifts in certain situations. To see why this is the case
# consider convolution by a boxcar kernel with an even width: Such kernels are
# intrinsically asymmetric... Recommend setting det_subgrid_fac to an odd
# number over 1.

import numpy as np
from EMToolKit.schemas import util
from EMToolKit.schemas.util import multivec_matmul, rindices, ftp, btp

class element_grid:
    # Set up the elements. Required inputs:
    #     grid: A coord_grid.
    #     parms: Parameters or anything else used to evaluate the response function.
    #     func: The response function evaluator.
    # Optional inputs:
    #     footprint: How far away from the input point to evaluate the response functions,
    #                in grid points. Default: at least 21 points.
    #     stencil_thold: In addition to the footprint, a stencil is computed to
    #                    determine which grid points to use evaluate around the input point.
    #                    An initial check of the response function (with evaluation point
    #                    at the origin) is made, and points falling below this threshold are
    #                    omitted. Would probably be better to check the stencil for every
    #                    element's location, although that would be slower... Default: .0005
    #     nsubgrid: To take into accound subgrid-scale effects, evaluate the response/basis functions.
    #               at this multiple of the grid scale. Default: 3
    #     thold: Final threshold for keeping points when evaluating the response/basis functions.
    #            Default: 0.005.
    def __init__(self, grid, parms, func, footprint=None, stencil_thold=5.0e-4, nsubgrid=3, thold=0.005):
        [self.thold, self.coords, self.parms, self.func] = [thold, grid, parms, func]
        [self.nelm, self.nsubgrid] = [np.prod(grid.dims), np.broadcast_to(nsubgrid,grid.dims.shape)]
        [self.subgrid, self.eval_grid] = [grid.subgrid(fac=self.nsubgrid), self.get_eval_grid()]
        # Generate the stencil:
        if(footprint is None): fpoffset = np.ceil(10*self.nsubgrid/2).astype(np.int32)
        else: fpoffset = np.ceil((footprint/self.nsubgrid-self.nsubgrid)/2).astype(np.int32)
        self.stencil = (rindices(self.nsubgrid+2*fpoffset) - fpoffset - 0.5*(nsubgrid-1.0))
        vals = self.evaluate(self.coords.origin)[0].flatten()
        self.stencil = np.vstack([x.flatten()[vals >= stencil_thold] for x in list(ftp(self.stencil))]).T
        # Set the number of addresses:
        self.nadr = self.get_nadr()
    # Standard grid for evaluating the response/basis functions is the same as the
    # subgrid. Note: indices for the eval_grid are assumed to be the same as the subgrid.
    def get_eval_grid(self): return self.coords.subgrid(fac=self.nsubgrid)
    # Standard assumption is number of addresses is same as number of elements:
    def get_nadr(self): return self.nelm
    # Evaluate the source/basis function at a given point
    def evaluate(self, point):
        subpt = self.subgrid.inds(point) # Find where the point is relative to the subgrid
        # Find the stencil evaluation indices (which are registered to the subgrid)
        # in the vicinity of this point.
        subinds = np.round(self.stencil+subpt+tiny)
        # Get the coordinates of these evaluation indices in the evaluation coordinate frame
        # and the subgrid coordinates:
        [eval_coords, output_coords] = [self.eval_grid.coords(subinds), self.subgrid.coords(subinds)]
        # Compute the response of these evaluation points to the input point, and their
        # coordinates:
        return self.func(self.eval_grid.coords(subpt), eval_coords, self.parms), output_coords
    # Return the elements addressed by a given index:
    def elements(self,index):
        # Map the address to a coordinate and run the element_grid's evaluator:
        [vals,coords] = self.evaluate(self.coords.coords(np.array(np.unravel_index(index,self.coords.dims))))
        # Return the result, element index is same as input address:
        return index+0*vals.astype(np.int32), vals, coords
    # Run the evaluator for the given point and compute the output value and indices to flatinds:
    def response(self,point): return self.coords.flatinds(*self.evaluate(point), thold=self.thold)

# The detector grid is a straight implementation of the base class:
class detector_grid(element_grid): pass

# The source grid is the same as the base class except that the evaluation grid for
# the source basis functions is a subgrid consisting of indices rather than using a
# physically dimensioned coordinate system. We use a fun trick with coord_grid's
# identify and subgrid methods to create this.
class source_grid(element_grid):
    def get_eval_grid(self): return self.coords.identity().subgrid(fac=self.nsubgrid)

