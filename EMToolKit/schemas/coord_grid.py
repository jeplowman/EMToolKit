# There is an issue with rounding and floating point jitter for aligned grids that
# share some ratios in their spacing. We add this small offset when discretizing
# grid indices to avoid this, which is a bit of a hack...
tiny = 1.0e-4

import numpy as np
from EMToolKit.schemas.util import multivec_matmul, rindices, ftp, btp

# This implements the notion of a grid of coordinate points that are meant to be indexed by integers.
# The basic class includes implementation for an affine (linear transformation) grid, but more general
# coordinate systems will not be hard to extend. The transformation from the indices into the coordinate
# points here should be reversible or the inds and flatinds methods will not work.
class coord_grid:
    # Set up the grid. Inputs:
    #    dims: The dimensions of the source grid, a 1-D array containing [nx0,nx1,...]
    #    origin: The location of the center of the [0,0,...] pixel in the grid
    #    fwd: The information needed to define the forward transformation (i.e., from
    #         indices to coordinates), in addition to the origin. For the affine
    #         transformation, this is an nd by nd matrix, where nd is the number of dimensions.
    #    inv: The information needed to define the inverse transformation (i.e., from
    #         coordinates to indices). Optional -- For the base (affine) implementation, this is
    #         computed from fwd using np.linalg.inv.
    #    frame: Information about the coordinate frame. Not used by the coord_grid itself, but
    #           may be used by other routines that want to transform between grids.
    def __init__(self, dims, origin, fwd, frame, inv=None):
        [self.dims,self.origin,self.fwd,self.inv,self.frame] = [dims,origin,fwd,inv,frame]
        if(self.inv is None): self.inv = self.get_inv()
    # Routine to set up the parameters of the grid inverse (from coordinates to indices)
    # operation. This default assumes an affine transformation and uses the matrix inverse.
    def get_inv(self): return np.linalg.inv(self.fwd)
    # Get a grid that's clocked to this grid, but is higher resolution by an integer factor.
    def subgrid(self,fac=2):
        return coord_grid(self.dims*fac,self.coords(0.5/fac-0.5+0.0*self.dims),self.fwd/fac,self.frame)
    # Get the identify grid for this grid -- i.e., it maps indices to themselves. Indices are coordinates
    # too, ya know! Pairs well with subgrid.
    def identity(self): return coord_grid(self.dims,0.0*self.dims,np.diag(1+0.0*self.dims),np.arange(len(self.dims)))
    # Returns indices given a set of coordinates. Does no discretize for various reasons.
    # Order is reversed, and the inv operator transposed, due to how numpy array broadcasting
    # and matrix operations interact. Because coords returns coordinates reference to the
    # centers of the grid elements, discretization of these indices should be done with
    # rounding, not flooring (see flatinds). The domain of each grid element extends
    # 0.5 grid spacings to each side.
    def inds(self,coords): return multivec_matmul(self.inv,coords-self.origin)
    # Returns coordinates given a set of indices. Coordinates returned for an integer index
    # are for the center of the grid element, not its corner.
    def coords(self,inds): return multivec_matmul(self.fwd, inds)+self.origin
    # Returns the 'flattened' indices given a set of coordinates. Does discretize (because it has to).
    # Also discards out-of-bounds points. Because of this, there's an accompanying vals array that can
    # be used to account for the discarding.
    def flatinds(self,vals,coords,thold=0):
        [inds, keeps] = [list(np.round(self.inds(coords)+tiny).T.astype(np.int32)), vals>thold]
        for j in range(0,len(self.dims)): keeps *= (inds[j] >= 0)*(inds[j] < self.dims[j])
        return vals[keeps], np.ravel_multi_index(inds,self.dims,mode='clip')[keeps]
