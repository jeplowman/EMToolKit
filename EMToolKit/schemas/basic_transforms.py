# This is a minimal implementation of the coordinate transform used by get_sparse_response_matrix.
# All it does is compare the 'names' attributes of the input coords to that of the output coords.
# This works for identical coordinates, rearrangements, and downprojections, but anything
# more complex will get zeroed out. In terms of the necessary pieces for instantiating
# and applying it, though, everything should work the same from the outside and it can be
# subclassed to add the features for any desired transformation. Note that
# these coordinate transforms might NOT be invertible in general (e.g., projection)
# unlike the transforms in the coord_grid class.

import numpy as np
from EMToolKit.schemas.util import multivec_matmul

class generic_transform:
    # Setup: takes two objects defining the two coordinate system. The minimal implementation
    # just compares the names attributes. This can be changed by overriding the
    # init_transform() method. Similarly, the transform itself is taken to be
    # an affine transform, but this can be changed by overriding the transform method.
    def __init__(self,coords_in,coords_out):
        [self.coords_in,self.coords_out] = [coords_in, coords_out]
        self.init_transform()
    def init_transform(self):
        [ndin,ndout] = [len(self.coords_in.frame.names), len(self.coords_out.frame.names)]
        [self.origin,self.fwd] = [np.zeros(ndin), np.zeros([ndout,ndin])]
        for i in range(0,ndout):
            for j in range(0,ndin): self.fwd[i,j] = self.coords_out.frame.names[i] == self.coords_in.frame.names[j]
    # As currently written will only work on one coordinate point at at time...
    def transform(self,coords): return multivec_matmul(self.fwd,coords)+self.origin

# A trivial coordinate frame object for use by the minimal implementation of the coord_transform
# class. Only contents are a set of names for the coordinates.
class trivialframe:
    def __init__(self,names): self.names=names

# This is a very basic transform between fits array index and physical
# (e.g, small-angle plane-of-sky) coordinate systems. It's here for a
# demo and placeholder for future implementation using the astropy coordinates object.
# That future implementation should behave equivalently to this when used on simple
# examples.
class fits_transform(object):

    def __init__(self,header):
        self.crpix1, self.crpix2 = header.get('CRPIX1',1), header.get('CRPIX2',1)
        self.crval1, self.crval2 = header.get('CRVAL1',0), header.get('CRVAL2',0)
        self.cdelt1, self.cdelt2 = header.get('CDELT1',1), header.get('CDELT2',1)
        self.crota = header.get('CROTA',None)
        if(self.crota) is None: self.crota = header.get('CROTA2',0.0)

        theta = self.crota*np.pi/180.0
        self.tfmat = np.zeros([2,2])
        self.tfmat[0,0] = np.cos(theta)
        self.tfmat[0,1] = -np.sin(theta)
        self.tfmat[1,0] = np.sin(theta)
        self.tfmat[1,1] = np.cos(theta)

    def index2coord(self,i,j):

        x0 = (i-(self.crpix1-1))*self.cdelt1
        y0 = (j-(self.crpix2-1))*self.cdelt2
        theta = self.crota*np.pi/180.0
        x = x0*np.cos(theta) - y0*np.sin(theta)
        y = x0*np.sin(theta) + y0*np.cos(theta)
        return np.array([x+self.crval1, y+self.crval2])

    def coord2index(self,xr,yr):

        theta = self.crota*np.pi/180.0
        x = xr-self.crval1
        y = yr-self.crval2
        x0 = x*np.cos(-theta) - y*np.sin(-theta)
        y0 = x*np.sin(-theta) + y*np.cos(theta)
        return np.array([x0/self.cdelt1+(self.crpix1-1), y0/self.cdelt2+(self.crpix2-1)])
