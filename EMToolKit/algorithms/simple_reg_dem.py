#+
# Just repeating the IDL comment block here for now:
# [dems,chi2] = simple_reg_dem(data, errors, exptimes, logt, tresps,
#					kmax=100, kcon=5, steps=[0.1,0.5], drv_con=8.0, chi2_th=1.0, tol=0.1):
#
# Computes a DEM given a set of input data, errors, exposure times, temperatures, and
# temperature response functions. Inputs:
#	data: image data (dimensions n_x by n_y by n_channels)
#	errors: image of uncertainties in data (dimensions n_x by n_y by n_channels)
#	exptimes: exposure times for each image (dimensions n_channels)
#	logt: Array of temperatures (dimensions n_temperatures)
#	tresps: Arrays of temperature response functions (dimensions n_temperatures by n_channels)
#
# Returns:
#	array of DEMs (dimensions n_x by n_y by n_temperatures). Dimensions depend on the units
#	of logt (need not be log temperature, despite the name) and tresps. For instance,
#	if tresps has units cm^5 and the logt values are base 10 log of temperature in Kelvin,
#	the output DEM will be in units of cm^{-5} per unit Log_{10}(T).
#
# Outputs:
#	chi2: Array of reduced chi squared values (dimensions n_x by n_y)
#
# Optional Keywords:
#	kmax: The maximum number of iteration steps (default - 100).
#	kcon: Initial number of steps before terminating if chi^2 never improves (default - 5).
#	steps: Two element array containing large and small step sizes (default - 0.1 and 0.5).
#	drv_con: The size of the derivative constraint - threshold delta_0 limiting the change in
#			 log DEM per unit input temperature (default - 8; e.g., per unit log_{10}(T)).
#	chi2_th: Reduced chi^2 threshold for termination (default - 1).
#	tol: How close to chi2_th the reduced chi^2 needs to be (default - 0.1).
#
# The iteration for each pixel is compiled with Numba and the pixels are run in parallel
# on every core. It reproduces the original NumPy/SciPy loop to rounding error (within
# 4e-10 relative on the AIA example data, including off-limb pixels with no signal).
# The first call compiles it, which takes several seconds, and caches it to disk for
# later sessions.
#
# Author: Joseph Plowman -- 09-15-2021
# See: https://ui.adsabs.harvard.edu/abs/2020ApJ...905...17P/abstract
#-
import numba
import numpy as np


def simple_reg_dem(data: np.ndarray, errors: np.ndarray, exptimes: np.ndarray, logt: np.ndarray, tresps: np.ndarray,
					kmax: int = 100, kcon: int = 5, steps: list[float] = [0.1,0.5], drv_con: float = 8.0,
					chi2_th: float = 1.0, tol: float = 0.1) -> tuple[np.ndarray, np.ndarray]:

	[nt,nd] = tresps.shape
	nt_ones = np.ones(nt)
	[nx, ny, nd] = data.shape
	dT = logt[1:nt]-logt[0:nt-1]
	[dTleft, dTright] = [np.diag(np.hstack([dT,0])),np.diag(np.hstack([0,dT]))]
	[idTleft, idTright] = [np.diag(np.hstack([1.0/dT,0])),np.diag(np.hstack([0,1.0/dT]))]
	Bij = ((dTleft+dTright)*2.0 + np.roll(dTright,-1,axis=0) + np.roll(dTleft,1,axis=0))/6.0
	Rij = np.matmul((tresps*np.outer(nt_ones,exptimes)).T,Bij) # Matrix mapping coefficents to data
	Dij = idTleft+idTright - np.roll(idTright,-1,axis=0) - np.roll(idTleft,1,axis=0)
	regmat = Dij*nd/(drv_con**2*(logt[nt-1]-logt[0]))
	rvec = np.sum(Rij,axis=1)

	npix = nx*ny
	dems = np.zeros([npix,nt])
	chi2 = np.zeros(npix) - 1.0
	_simple_reg_dem_pixels(
		np.ascontiguousarray(data, dtype=np.float64).reshape(npix, nd),
		np.ascontiguousarray(errors, dtype=np.float64).reshape(npix, nd),
		np.ascontiguousarray(Rij, dtype=np.float64), np.ascontiguousarray(regmat, dtype=np.float64),
		np.ascontiguousarray(rvec, dtype=np.float64), int(kmax), int(kcon), float(steps[0]), float(steps[1]),
		float(chi2_th), float(tol), max(1, min(npix, numba.get_num_threads())), dems, chi2)

	return dems.reshape([nx,ny,nt]), chi2.reshape([nx,ny])


# error_model="numpy" makes a division by zero give an infinity or NaN instead of raising,
# as NumPy does. "contract" and "reassoc" let the compiler fuse and reorder the sums, which
# is about a third faster and changes the result by ~1e-12 relative. The other fast-math
# flags would let the compiler drop the NaN checks that decide when a pixel has failed.
_flags = dict(error_model="numpy", fastmath={"contract", "reassoc"})


@numba.njit(**_flags)
def _chi2(dat0: np.ndarray, ierr: np.ndarray, Rij: np.ndarray, s: np.ndarray, exp_s: np.ndarray) -> float:
	"""Reduced chi^2 of the coefficients exp(s), leaving exp(s) in exp_s."""
	nd, nt = Rij.shape
	for j in range(nt): exp_s[j] = np.exp(s[j])
	total = 0.0
	for i in range(nd):
		model = 0.0
		for j in range(nt): model += Rij[i,j]*exp_s[j]
		residual = (dat0[i]-model)*ierr[i]
		total += residual*residual
	return total/nd


@numba.njit(**_flags)
def _cho_factor(a: np.ndarray) -> bool:
	"""
	Overwrite the lower triangle of a with its Cholesky factor. False where
	scipy.linalg.cho_factor would raise: a is not finite or not positive definite.
	"""
	n = a.shape[0]
	for j in range(n):
		d = a[j,j]
		for k in range(j): d -= a[j,k]*a[j,k]
		if not (d > 0.0) or not np.isfinite(d): return False
		d = np.sqrt(d)
		a[j,j] = d
		for i in range(j+1, n):
			x = a[i,j]
			for k in range(j): x -= a[i,k]*a[j,k]
			a[i,j] = x/d
	return True


@numba.njit(**_flags)
def _cho_solve(c: np.ndarray, b: np.ndarray, x: np.ndarray) -> None:
	"""Solve A x = b, given the lower Cholesky factor c of A."""
	n = c.shape[0]
	for i in range(n):
		v = b[i]
		for k in range(i): v -= c[i,k]*x[k]
		x[i] = v/c[i,i]
	for i in range(n-1, -1, -1):
		v = x[i]
		for k in range(i+1, n): v -= c[k,i]*x[k]
		x[i] = v/c[i,i]


@numba.njit(**_flags)
def _simple_reg_dem_pixel(data: np.ndarray, err: np.ndarray, Rij: np.ndarray, regmat: np.ndarray, rvec: np.ndarray,
							kmax: int, kcon: int, step0: float, step1: float, chi2_th: float, tol: float,
							dem: np.ndarray, s: np.ndarray, exp_s: np.ndarray, s_trial: np.ndarray, exp_trial: np.ndarray,
							deltas: np.ndarray, rhs: np.ndarray, solution: np.ndarray, gram: np.ndarray, amat: np.ndarray,
							dat0: np.ndarray, ierr: np.ndarray, dat: np.ndarray) -> float:
	"""
	The loop body of the original simple_reg_dem for one pixel. Writes the DEM into dem and
	returns chi^2. The arguments after dem are work space, overwritten here.
	"""
	nd, nt = Rij.shape

	# A zero or NaN error gives a NaN initial guess in the original, so cho_factor raises on
	# the first step, and the pixel is left with a NaN DEM and chi^2 of -1.
	for i in range(nd):
		if err[i] == 0.0 or err[i] != err[i]:
			for j in range(nt): dem[j] = np.nan
			return -1.0

	numerator = 0.0
	denominator = 0.0
	for i in range(nd):
		x = data[i]
		dat0[i] = x if (x != x or x > 0.0) else 0.0 # np.clip(data, 0.0, None), which keeps NaN
		ierr[i] = 1.0/err[i]
		floored = dat0[i] if (dat0[i] != dat0[i] or dat0[i] > 1.0e-2) else 1.0e-2
		numerator += rvec[i]*(floored/err[i]**2)
		denominator += (rvec[i]/err[i])**2
	ratio = numerator/denominator
	# With "reassoc", exp(log(x)) compiles to x, which would turn the NaN of the log of a
	# negative number back into a number to iterate on.
	if not ratio >= 0.0: ratio = np.nan
	s0 = np.log(ratio)
	for j in range(nt):
		s[j] = s0
		exp_s[j] = np.exp(s0)

	# (Rij/err)^T (Rij/err), which is the same for every step
	for j in range(nt):
		for k in range(j+1):
			x = 0.0
			for i in range(nd): x += Rij[i,j]*ierr[i]*(Rij[i,k]*ierr[i])
			gram[j,k] = x

	chi2 = -1.0
	for k in range(kmax):
		# amat = mmat^T mmat + regmat, lower triangle only
		for j in range(nt):
			for l in range(j+1): amat[j,l] = exp_s[j]*gram[j,l]*exp_s[l] + regmat[j,l]
		if not _cho_factor(amat): break

		c2p = 0.0
		for i in range(nd):
			model = 0.0
			model_linear = 0.0
			for j in range(nt):
				model += Rij[i,j]*exp_s[j]
				model_linear += Rij[i,j]*((1.0-s[j])*exp_s[j])
			residual = (dat0[i]-model)*ierr[i]
			c2p += residual*residual
			dat[i] = (dat0[i]-model_linear)/err[i] # Correct data by f(s)-s*f'(s)...
		c2p /= nd

		for j in range(nt):
			x = 0.0
			for i in range(nd): x += Rij[i,j]*ierr[i]*dat[i]
			rhs[j] = x*exp_s[j]
		_cho_solve(amat, rhs, solution)

		# A NumPy float, so that with Numba disabled a zero maximum gives NaN as it does compiled
		largest = np.float64(0.0)
		for j in range(nt):
			deltas[j] = solution[j]-s[j]
			d = abs(deltas[j])
			if d > largest or d != d: largest = d
		limit = 0.5/step0
		scale = (largest if largest < limit or largest != largest else limit)/largest
		for j in range(nt): deltas[j] *= scale

		ds = -1.0 if c2p < chi2_th else 1.0 # Direction sign; is chi squared too large or too small?
		for j in range(nt): s_trial[j] = s[j]+deltas[j]*ds*step0
		c20 = _chi2(dat0, ierr, Rij, s_trial, exp_trial)
		for j in range(nt): s_trial[j] = s[j]+deltas[j]*ds*step1
		c21 = _chi2(dat0, ierr, Rij, s_trial, exp_trial)
		interp_step = (step0*(c21-chi2_th)+step1*(chi2_th-c20))/(c21-c20)
		if interp_step < step0: interp_step = step0 # np.clip, which keeps NaN
		elif interp_step > step1: interp_step = step1
		for j in range(nt): s[j] += deltas[j]*ds*interp_step
		chi2 = _chi2(dat0, ierr, Rij, s, exp_s)

		if (ds*(c2p-c20)/step0 < tol and k > kcon) or abs(chi2-chi2_th) < tol: break

	for j in range(nt): dem[j] = np.exp(s[j])
	return chi2


@numba.njit(parallel=True, cache=True, **_flags)
def _simple_reg_dem_pixels(data: np.ndarray, errors: np.ndarray, Rij: np.ndarray, regmat: np.ndarray,
							rvec: np.ndarray, kmax: int, kcon: int, step0: float, step1: float, chi2_th: float,
							tol: float, nchunks: int, dems: np.ndarray, chi2: np.ndarray) -> None:
	"""
	Every pixel of data and errors, (n_pixels, n_channels), into dems, (n_pixels, n_temperatures),
	and chi2, (n_pixels,).

	The pixels are dealt into nchunks chunks, every nchunks-th pixel to the same chunk, so that
	each thread draws pixels from the whole image rather than one block of it (an active region
	in one block would otherwise keep one thread busy while the others sit idle). Each chunk
	allocates its work arrays once. nchunks is passed in, since calling numba.get_num_threads()
	here would stop Numba from caching this function.
	"""
	npix, nd = data.shape
	nt = Rij.shape[1]
	for c in numba.prange(nchunks):
		s = np.empty(nt)
		exp_s = np.empty(nt)
		s_trial = np.empty(nt)
		exp_trial = np.empty(nt)
		deltas = np.empty(nt)
		rhs = np.empty(nt)
		solution = np.empty(nt)
		gram = np.empty((nt, nt))
		amat = np.empty((nt, nt))
		dat0 = np.empty(nd)
		ierr = np.empty(nd)
		dat = np.empty(nd)
		for p in range(c, npix, nchunks):
			chi2[p] = _simple_reg_dem_pixel(data[p], errors[p], Rij, regmat, rvec, kmax, kcon, step0, step1,
											chi2_th, tol, dems[p], s, exp_s, s_trial, exp_trial, deltas, rhs,
											solution, gram, amat, dat0, ierr, dat)
