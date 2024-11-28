import time
import resource
import numpy as np
#from processify import processify

from scipy.sparse.linalg import LinearOperator

# These objects compute the forward operator for a combined spatial and temperature
# (DEM) forward problem, along with it's transpose. Given a proposed solution
# vector (e.g., an [T, x, y] cube) in both temperature and space (dimensioned 
# [ntemp, nx, ny] or [ntemp, nx*ny]) the forward operator computes the data vector,
# which has arranged as a flat vector [d0_0, d0_1, ... d0_n0, d1_0, d1_1, ... d1_n1,
# ... dm_0, dm_1, ... dm_nm], where the data is organized into m 'channels', all
# with the same temperature response function (e.g., the EUV channels in AIA). The
# number of data elements in each channel (n0, ..., nm) need not all be the same.
#
# The transpose of the forward operator acts on a data vector; it gives the sum of
# all of the spatio-thermal response functions in the forward operator, weighted by
# their values in the data vector. For example, a pixel in AIA 171 has a response
# function that consists of a spatial PSF (plus pixelization) which responds spatially to
# a small neighborhood of the source cube and a temperature response function
# which responds thermally according to the 171 temperature response functions.
# If the transpose operates on a data vector with only that pixel and channel 'lit up',
# The output is a cube in the source space with the spatial response function of
# that pixel lit according to the PSF at the temperatures 171 responds to. All other
# source cube elements will be zero.
#
# The outputs of the forward and transpose operations are flattened since that's what
# the solvers expect.
class multi_instrument_linear_operator(LinearOperator):
	def __init__(self, operators, wrapargs={}):
		self.dtype=None
		self.rm_timer = 0
		self.lm_timer = 0
		self.nchans = len(operators)
		self.operators = operators
		self.outrange_lo = np.zeros(self.nchans,dtype=np.uint64)
		self.outrange_hi = np.zeros(self.nchans,dtype=np.uint64)
		self.nb_temp = operators[0].nb_temp
		self.nb_spat = operators[0].nb_spat
		self.norms = wrapargs.get('norms',np.ones(self.nchans))
		for i in range(0,self.nchans):
			if(i > 0): self.outrange_lo[i] = self.outrange_hi[i-1]
			self.outrange_hi[i] = self.outrange_lo[i]+operators[i].shape[0]
			if(operators[i].nb_temp != self.nb_temp or operators[i].nb_spat != self.nb_spat):
				print('Warning: multi_instrument_linear_operator inputs have differing dimensions.')
		self.ndat = np.max(self.outrange_hi)
		self.nsrc = operators[0].shape[1]
		self.shape = [self.ndat,self.nsrc]
		
	def _matvec(self, rightvec_in):
		tstart = time.time()
		outvec = np.zeros(self.ndat)
		for i in range(0,self.nchans):
			outvec[self.outrange_lo[i]:self.outrange_hi[i]] = (self.operators[i]*rightvec_in)/self.norms[i]
		self.rm_timer += time.time()-tstart
		return outvec
		
	def _rmatvec(self, leftvec_in):
		tstart = time.time()
		outvec = np.zeros(self.nsrc)
		for i in range(0,self.nchans):
			outvec += (self.operators[i].T*leftvec_in[self.outrange_lo[i]:self.outrange_hi[i]])/self.norms[i]
		self.lm_timer += time.time()-tstart
		return outvec


# This represents a forward operator where the spatial and thermal response
# separable. Most imagers with a fixed EUV channel are like this, although
# the diffraction spikes in AIA break it slightly. If the temperatures
# are input (in the 'temps' keyword), this version assumes that tresp is a 
# conventional tabulated temperature response function and that 'temps' lists
# the corresponding temperatures as inputs. It then assumes the temperature
# basis elements are 'spike' functions as in Plowman & Caspi 2020. This is needed
# to be consistent with the derivative-based regularization operators found
# elsewhere in this code. If the 'temps' keyword is not supplied, tresp is 
# assumed to have the basis elements already rolled into it (supplying 
# a raw temperature response function this way will be like a top hat 'bin'
# function but will be missing the bin width unless that is premultipied).
# The spatial response matrix must be input separately; it can be computed
# using the element function framework used for the SPICE PSF correction.
class single_instrument_linear_operator_separable(LinearOperator):
	def __init__(self, sresp, tresp, temps=None, exptime=None):
		self.dtype=None
		self.nb_spat = sresp.shape[1]
		self.nb_temp, self.nchans = tresp.size, 1
		self.sresp = sresp
		self.tresp = tresp
		self.nsrc = self.nb_spat*self.nb_temp
		self.shape = [sresp.shape[0],self.nsrc]
		if(tresp.ndim==1): tresp = np.expand_dims(tresp,1)
		if(temps is not None):
			nt = len(temps)
			nt_ones = np.ones(nt)
			dT = temps[1:nt]-temps[0:nt-1]
			if(exptime is None): exptime = 1.0
			if(np.isscalar(exptime)): exptimes = np.array([exptime])
			else: exptimes = exptime
			[dTleft, dTright] = [np.diag(np.hstack([dT,0])),np.diag(np.hstack([0,dT]))]
			Bij = ((dTleft+dTright)*2.0 + np.roll(dTright,-1,axis=0) + np.roll(dTleft,1,axis=0))/6.0
			# This is from the multidimensional version. For a single temperature response
			# it can be simplified; left here because it's not important enough to spool up
			# and bugfix right now.
			self.Rij = np.matmul((tresp*np.outer(nt_ones,exptimes)).T,Bij) # Matrix mapping coefficents to data
		else:
			self.Rij = tresp
		self.Rij = self.Rij.flatten()
		# This makes Rij the right shape for numpy to expand it
		# when multiplied by another vector on the left.
		self.Rij_expanded = np.expand_dims(self.Rij,1) 

	# Multiply on the right of the operator (the conventional forward way):
	def _matvec(self, rightvec_in):
		rightvec = rightvec_in.reshape([self.nb_temp,self.nb_spat])		
		outvec = np.zeros(self.nb_spat)
		for j in range(0,self.nb_temp):
			outvec += self.Rij[j]*rightvec[j]
		return self.sresp*outvec

	# Multiply a vector on the left of the operator (the transpose): 
	def _rmatvec(self, leftvec_in):
		return ((self.sresp.T*leftvec_in)*self.Rij_expanded[:,:]).flatten()


# The regularization operators shown below do not need a distinct
# left multiplication method since the regularization is assumed to be
# symmetric. Regularizations map from the source coefficient space
# back to the source coefficient space, whereas the forward operator
# maps from source coefficient space to data number space. The ATA matrix
# that's part of the chi squared formulation turns the forward operator
# into an operator that maps source coefficients to source coefficients --
# it multiplies the transpose of the forward operator by the forward
# operator, which is why the transpose is needed. The ATA matrix is not
# explicitly computed, however, since it is much larger than its
# component parts.

## This wraps a regularization matrix for use by the sparse nlmap solver.
#class reg_operator_wrapper(LinearOperator):
#	def __init__(self,regmat):
#		self.regmat=regmat
#		self.shape=regmat.shape
#		self.dtype=regmat.dtype
#
#	def _matvec(self,vec):
#		outvec = self.regmat*vec
#		return outvec

# This wrapper is for sparse_d2_partial_matrix with post-factor:		
class reg_operator_postfac_wrapper(LinearOperator):
	def __init__(self,operator):
		self.operator = operator
		self.shape = self.operator[0].shape
		self.dtype = None

	def _matvec(self,vec):
		return self.operator[1]*self.operator[0]*vec
        
	def _rmatvec(self,vec):
		return self.operator[1]*self.operator[0].T*vec

# This computes the regularization matrix for minimizing the derivative along
# a spatial direction, like the temperature derivative matrix in Plowman & Caspi 2020.
# This version uses a diagonal sparse matrix format which should be more space efficient.
# The use_postfactor argument allows the matrix to be split into a post factor part and
# the actual matrix. This can lead to further space savings if the values of the matrix
# are all the same. However it returns the two pieces rather than the single matrix
# so the outputs need to be treated differently. Changing this behavior would require
# implementing the output as a full sparse matrix object including the glue to
# work with numpy/scipy's operator overloading.
def sparse_d2_partial_matrix(dims, axis, nd, steps=None, drv_con=8.0, dtype=np.float32, use_postfactor=False):
	from scipy.sparse import diags, csr_matrix

	if(use_postfactor==True): dtype_internal = np.int8
	if(use_postfactor==False): dtype_internal = dtype
    
	npt = np.prod(dims)
	if(steps is None): steps = np.ones(dims[axis]-1)

	# Easiest way to update all elements by a specified axis
	# which is not know ahead of time is to transpose it so the axis comes
	# first, update that way, then transpose back.
	tparg = np.arange(len(dims),dtype=np.int32)
	tparg[axis] = tparg[0]
	tparg[0] = axis

	# Diagonal terms:
	diagvec = np.zeros(dims,dtype=dtype_internal).transpose(tparg)
	for i in range(0,dims[axis]-1):
		if(use_postfactor):
			diagvec[i] += 1
			diagvec[i+1] += 1
		else:
			diagvec[i] += 1.0/steps[i]
			diagvec[i+1] += 1.0/steps[i]
	diagvec = diagvec.transpose(tparg)
	
	diagvec = diagvec.flatten()

	offdiagvec = np.zeros(dims).transpose(tparg).T
	if(use_postfactor==False): offdiagvec[...,0:-1] = -1.0/steps
	if(use_postfactor==True): offdiagvec[...,0:-1] = -1
	offdiagvec = offdiagvec.T.transpose(tparg).flatten()
    
	if(axis < len(dims)-1): offset = np.prod(dims[axis+1:]).astype(np.int64)
	else: offset = 1

	output_matrix = diags((diagvec,offdiagvec,offdiagvec), offsets = (0,-offset,offset), shape=(npt,npt), dtype=dtype_internal, format='dia')
    
	if(use_postfactor):
		output = output_matrix,(nd/(np.sum(steps)*drv_con**2)/np.median(steps)).astype(dtype)
		return reg_operator_postfac_wrapper(output)
	else:
		return output_matrix*nd/(np.sum(steps)*drv_con**2)

