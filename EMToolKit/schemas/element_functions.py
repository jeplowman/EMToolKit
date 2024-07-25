import numpy as np
from EMToolKit.schemas.util import multivec_matmul
tiny = 1.0e-4

# This is a basis function that's a set of rectilinear 'in' or 'out' boxes -- e.g., pixels in 2D:
def bin_function(ptcoords,coordarr,parms): return np.prod(np.round((ptcoords-coordarr).T+tiny) == 0,axis=0)

# Rectilinear spike function. Can be 1-D or multidimensional. Should work well for a
# linear (or bilinear for 2D, etc) interpolating basis function. Note: Test this...
def spike_function(ptcoords,coordarr,parms): return np.prod(np.clip(1.0-np.abs(ptcoords-coordarr).T,0,1),axis=0)

# This is a setup function for a 3D or 2D Gaussian type PSF/response function.
# It's defined in terms of the 3 axes of the ellipse and a set of three angles
# about which the ellipse is rotated. The initial axes are x, y, z, while
# the rotation matrices use the Tait Bryan convention z, x, y -- i.e., if
# the angles are all zero, sigmas[0] will be the ellipse length along x,
# sigmas[1] wil be the ellipse length along y, etc; and, if the angles are not zero
# the ellipse is first rotated about the z axis (in the x-y plane), then x (z-y plane)
# then y (z-x plane). If the PSF is evaluated in 2D, only the first two
# axis lengths and the first angle will be used.
def get_3d_cov(sigmas,angles): # Get covariance for 3D Mahalanobis distance:
    [[c1,c2,c3],[s1,s2,s3]] = [np.cos(angles),np.sin(angles)]
    # Compute covariance using Tait-Bryan angles ordered z, x, y:
    vec0 = sigmas[0]*np.array([c1*c3-s1*s2*s3, c3*s1+c1*s2*s3, -c2*s3])
    vec1 = sigmas[1]*np.array([-c2*s1, c1*c2, s2])
    vec2 = sigmas[2]*np.array([c1*s3+c3*s1*s2, s1*s3-c1*c3*s2, c2*c3])
    return np.outer(vec0,vec0)+np.outer(vec1,vec1)+np.outer(vec2,vec2)

def get_2d_cov(sigmas,theta):
	vec0 = sigmas[0]*np.array([np.cos(theta),np.sin(theta)])
	vec1 = sigmas[1]*np.array([-np.sin(theta),np.cos(theta)])
	return np.outer(vec0,vec0)+np.outer(vec1,vec1)

def spice_spectrograph_psf(pt, coords, inputs):
	slitwid = inputs[-1]
	if(len(inputs)<4): n_slit_subpts = 5
	else: n_slit_subpts = inputs[-2]
	subpts = np.zeros([n_slit_subpts,2])
	subpts[:,1] = slitwid*np.arange(n_slit_subpts)/n_slit_subpts - 0.5*slitwid + 0.5*slitwid/n_slit_subpts
	psf = nd_powgaussian_psf(pt+subpts[0], coords, inputs)
	for i in range(1,n_slit_subpts): psf += nd_powgaussian_psf(pt+subpts[i], coords, inputs)
	return psf/n_slit_subpts

# Evaluate an nd Gaussian PSF centered at the point pt, for each of the
# coordinates coords, based on the q (inverse of covariance) matrix q.
# Coords must have dimensions npts by nd where nd is the dimensionality
# of the Gaussian. q can be larger than nd by nd; higher dimensions will be ignored.
def nd_gaussian_psf(pt, coords, inputs):
    dxa = coords-pt
    q = inputs[0]
    return np.exp(-0.5*np.sum(dxa*multivec_matmul(q[0:pt.size,0:pt.size],dxa), axis=-1))

def nd_voigt_psf(pt, coords, inputs):
    from scipy.special import voigt_profile
    dxa = coords-pt
    q = inputs[0]
    g = inputs[1]*(np.log(2))**0.5
    e = inputs[2]
    mdist = np.sum(dxa*multivec_matmul(q[0:pt.size,0:pt.size],dxa), axis=-1)
    return voigt_profile(mdist**0.5,e,g)*(1.0/voigt_profile(0,e,g))


def nd_powgaussian_psf(pt, coords, inputs):
    dxa = coords-pt
    q = inputs[0]
    exp = inputs[1]
    return np.exp(-0.5*np.sum(dxa*multivec_matmul(q[0:pt.size,0:pt.size],dxa), axis=-1)**exp)

def flattop_guassian_psf(pt, coords, inputs):
	dxa = coords-pt
	q = inputs[0]
	flatness = inputs[1]
	flatsds = inputs[2]
	mdist = np.sum((dxa*multivec_matmul(q[0:pt.size,0:pt.size],dxa)), axis=-1)
	return (np.exp(-0.5*mdist)/(1.0-0.5*flatness*mdist*(np.exp(-0.5*mdist/(flatsds**2)))**2))
