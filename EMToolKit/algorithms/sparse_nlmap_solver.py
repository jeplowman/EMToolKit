import time
import resource
import numpy as np
from tqdm import tqdm

# Sparse in this usage here means that it uses sparse matrix methods,
# Not (necessarily) that it uses 'sparsity' of the solution as a constraint
# a la the L1 norm used in Mark Cheung's sparse_em basis pursuit algorithm.
# This is a result of a name collision in the mathematical terminology that
# we can't really avoid here. Sparse is the natural term for when a linear
# operator has mostly zero entries, and the accepted mathematical one.
# And although the solver here can work on non-sparse operators, its
# power is in very large dimensional problems which are only tractable
# when they're sparse.
from scipy.sparse.linalg import LinearOperator

# This operator implements the general linear operator for chi squared
# plus regularization with nonlinear mapping as outlined in Plowman &
# Caspi 2020.
class nlmap_operator(LinearOperator):
	def setup(self,amat,regmat,map_drvvec,wgtvec,reg_map_drvvec,dtype='float32',reg_fac=1):
		self.amat = amat
		self.regmat = regmat
		self.map_drvvec = map_drvvec
		self.wgtvec = wgtvec
		self.reg_map_drvvec = reg_map_drvvec
		self.dtype_internal=dtype
		self.reg_fac = reg_fac

	def _matvec(self,vec):
		chi2term = self.map_drvvec*(self.amat.T*(self.wgtvec*(self.amat*(self.map_drvvec*vec))))
		regterm = self.reg_map_drvvec*(self.reg_fac*self.regmat*(self.reg_map_drvvec*vec))
		return (chi2term+regterm).astype(self.dtype_internal)

	def _adjoint(self):
		return self

import time
import resource
import numpy as np

# Subroutine to do the inversion. Uses the nolinear mapping and iteration from Plowman & Caspi 2020 to ensure positivity of solutions.
# Can use nonlinear mappings other than logarithmic.
def solve(data0, errors0, amat0, guess=None, reg_fac=1, func=None, dfunc=None, ifunc=None, regmat=None, silent=False,
						solver=None, sqrmap=False, regfunc=None, dregfunc=None, iregfunc=None, map_reg=False, adapt_lam=True,
						solver_tol = 1.0e-3, niter=40, dtype='float32', steps=None, precompute_ata=False, flatguess=True, chi2_th=1.0,
						store_outer_Av=False, conv_chi2 = 1.0e-15):
	from scipy.sparse import diags
	from scipy.sparse.linalg import lgmres
	from scipy.linalg import cho_factor, cho_solve

	[ndat, nsrc] = amat0.shape

	# Being really careful that everything is the right dtype
	# so (for example) nothing gets promoted to double if dtype is single:
	solver_tol = np.dtype(dtype).type(solver_tol)
	zero = np.dtype(dtype).type(0.0)
	pt5 = np.dtype(dtype).type(0.5)
	two = np.dtype(dtype).type(2.0)
	one = np.dtype(dtype).type(1.0)
	conv_chi2 = np.dtype(dtype).type(conv_chi2)

	# A collection of example mapping functions.
	def idnfunc(s): return s # linear forward
	def iidnfunc(s): return s # linear inverse
	def didnfunc(s): return one + zero*s # linear derivative
	def expfunc(s): return np.exp(s) # exponential forward
	def iexpfunc(c): return np.log(c) # exponential inverse
	def dexpfunc(s): return np.exp(s) # exponential derivative
	def sqrfunc(s): return s*s # quadratic forward
	def isqrfunc(c): return c**pt5 # quadratic inverse
	def dsqrfunc(s): return two*s # quadratic derivative
	if(func is None or dfunc is None or ifunc is None):
		if(sqrmap): [func,dfunc,ifunc] = [sqrfunc,dsqrfunc,isqrfunc]
		else: [func,dfunc,ifunc] = [expfunc,dexpfunc,iexpfunc]
	if(regfunc is None or dregfunc is None or iregfunc is None):
		if(map_reg): [regfunc,dregfunc,iregfunc] = [idnfunc,didnfunc,iidnfunc]
		else: [regfunc,dregfunc,iregfunc] = [func,dfunc,ifunc]
	if(solver is None): solver = lgmres

	flatdat = data0.flatten().astype(dtype)
	flaterrs = errors0.flatten().astype(dtype)
	flaterrs[flaterrs == 0] = (0.05*np.nanmean(flaterrs[flaterrs > 0])).astype(dtype)

	guess0 = amat0.T*(np.clip(flatdat,np.min(flaterrs),None))
	guess0dat = amat0*(guess0)
	guess0norm = np.sum(flatdat*guess0dat/flaterrs**2)/np.sum((guess0dat/flaterrs)**2)
	guess0 *= guess0norm
	guess0 = np.clip(guess0,0.05*np.mean(np.abs(guess0)),None).astype(dtype)
	if(guess is None): guess = guess0
	if(flatguess): guess = ((1+np.zeros(nsrc))*np.mean(flatdat)/np.mean(amat0*(1+np.zeros(nsrc)))).astype(dtype)
	svec = ifunc(guess).astype(dtype)

	# This is an internal step length limiter to prevent overstepping
	# if the solution starts out in a highly nonlinear region of the mapping
	# function. The solver can still take larger steps since the limiter scales
	# it so that the smallest step size is maxdelta.
	maxdelta = iregfunc(np.max(guess0))-iregfunc(0.25*np.max(guess0))

	# Try these step sizes at each step of the iteration. Trial Steps are fast compared to computing
	# the matrix inverse, so having a significant number of them is not a problem.
	# Step sizes are specified as a fraction of the full distance to the solution found by the sparse
	# matrix solver (lgmres or bicgstab).
	if(steps is None): steps = np.array([0.00, 0.05, 0.15, 0.3, 0.5, 0.67, 0.85],dtype=dtype)
	minstep = np.min(steps[1:])
	nsteps = len(steps)
	step_loss = np.zeros(nsteps,dtype=dtype)

	reglam = one
	if(regmat is None and map_reg): regmat = diags(one/iregfunc(guess0)**two)
	if(adapt_lam and map_reg): reglam = (np.dot((regmat*svec),
							dfunc(svec)*(amat0.T*(1.0/flaterrs)))/
							np.dot((regmat*svec),(regmat*svec)))
	if(regmat is None and not map_reg):  regmat = diags(1.0/(guess0)**2)
	if(adapt_lam and not map_reg): reglam = (np.dot(dfunc(svec)*(regmat*guess),
							dfunc(svec)*(amat0.T*(1.0/flaterrs)))/
							np.dot(dfunc(svec)*(regmat*guess),dfunc(svec)*(regmat*guess)))

	# Still appears to be some issue with this regularization factor?
	# regmat = reg_fac*regmat*reglam # Old buggy application of reg_fac?
	reg_fac = reg_fac*reglam
	weights = (1.0/flaterrs**2).astype(dtype) # The weights are the errors...

	if(silent == False):
		print('Overall regularization factor:',reg_fac*reglam)

	if(not(precompute_ata)):
		nlmo = nlmap_operator(dtype=dtype,shape=(nsrc,nsrc))
		nlmo.setup(amat0,regmat,dfunc(svec),weights,dregfunc(svec),reg_fac=reg_fac)

	# --------------------- Now do the iteration:
	tstart = time.time()
	setup_timer = 0
	solver_timer = 0
	stepper_timer = 0
	for i in tqdm(range(0,niter), desc="Iteration"):
		tsetup = time.time()
		# Setup intermediate matrices for solution:
		dguess = dfunc(svec)
		dregguess = dregfunc(svec)
		bvec = dguess*amat0.T.dot(weights*(flatdat-amat0*(func(svec)-svec*dfunc(svec))))
		if(map_reg==False): bvec -= dregguess*(reg_fac*regmat*(regfunc(svec)-svec*dregfunc(svec)))
		setup_timer += time.time()-tsetup

		tsolver = time.time()
		# Run sparse matrix solver:
		[nlmo.map_drvvec,nlmo.reg_map_drvvec] = [dguess,dregguess]
		svec2 = solver(nlmo,bvec.astype(dtype),svec.astype(dtype),store_outer_Av=False,rtol=solver_tol.astype(dtype))
		svec2 = svec2[0]
		solver_timer += time.time()-tsolver

		tstepper = time.time()
		deltas = svec2-svec
		if(np.max(np.abs(deltas)) == 0): break # This also means we've converged.
		# Rescale the deltas so they don't exceed maxdelta at the smallest step size:
		deltas *= np.clip(np.max(np.abs(deltas)),None,maxdelta/minstep)/np.max(np.abs(deltas))

		# Try the step sizes:
		for j in range(0,nsteps):
			stepguess = func(svec+steps[j]*(deltas))
			stepguess_reg = regfunc(svec+steps[j]*(deltas))
			stepresid = (flatdat-amat0*(stepguess))*weights**pt5
			step_loss[j] = np.dot(stepresid,stepresid)/ndat + np.sum(stepguess_reg.T*(reg_fac*regmat*(stepguess_reg)))/ndat

		best_step = np.argmin((step_loss)[1:nsteps])+1 # First step is zero for comparison purposes...
		chi20 = np.sum(weights*(flatdat-amat0*(func(svec)))**two)/ndat
		reg0 = np.sum(regfunc(svec.T)*(reg_fac*regmat*(regfunc(svec))))/ndat

		# Update the solution with the step size that has the best Chi squared:
		svec = svec+steps[best_step]*(deltas)
		reg1 = np.sum(regfunc(svec.T)*(reg_fac*regmat*(regfunc(svec))))/ndat
		resids = weights*(flatdat-amat0*(func(svec)))**two
		chi21 = np.sum(weights*(flatdat-amat0*(func(svec)))**two)/ndat
		stepper_timer += time.time()-tstepper

		if(silent==False):
			print(round(time.time()-tstart,2),'s i =',i,'chi2 =',round(chi21,2),'step size =',round(steps[best_step],3), 'reg. param. =', round(reg1,2), 'chi2 change =',round(chi20-chi21,5), 'reg. change =',round(reg0-reg1,5))
			print('Setup: ', setup_timer, 'Solver: ', solver_timer, 'Stepper: ', stepper_timer)
			#print('lm_timer: ', amat0.lm_timer, 'rm_timer: ', amat0.rm_timer, 'reg_timer: ', regmat.rm_timer)
			print('New combined FOM:',chi21+reg1,'Old combined FOM:',chi20+reg0,'Change:',chi20+reg0-(chi21+reg1))
		if(np.abs(step_loss[0]-step_loss[best_step]) < conv_chi2 or chi21 < chi2_th): break # Finish the iteration if chi squared isn't changing

	return func(svec), chi21, resids
