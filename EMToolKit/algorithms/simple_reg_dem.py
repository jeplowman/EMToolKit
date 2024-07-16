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
# Author: Joseph Plowman -- 09-15-2021
# See: https://ui.adsabs.harvard.edu/abs/2020ApJ...905...17P/abstract
#-
def simple_reg_dem(data, errors, exptimes, logt, tresps,
					kmax=100, kcon=5, steps=[0.1,0.5], drv_con=8.0, chi2_th=1.0, tol=0.1):

	import numpy as np
	from scipy.linalg import cho_factor, cho_solve

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

	dems = np.zeros([nx,ny,nt])
	chi2 = np.zeros([nx,ny]) - 1.0
	for i in range(0,nx):
		for j in range(0,ny):
			err = errors[i,j,:]
			dat0 = np.clip(data[i,j,:],0.0,None)
			s = np.log(np.sum((rvec)*((np.clip(dat0, 1.0e-2, None))/err**2))/np.sum((rvec/err)**2)/nt_ones)
			for k in range(0,kmax):
				dat = (dat0-np.matmul(Rij,((1-s)*np.exp(s))))/err # Correct data by f(s)-s*f'(s)...
				mmat = Rij*np.outer(1.0/err,np.exp(s))
				amat = np.matmul(mmat.T,mmat)+regmat
				try: [c,low] = cho_factor(amat)
				except: break
				c2p = np.mean((dat0-np.dot(Rij,np.exp(s)))**2/err**2)
				deltas = cho_solve((c,low),np.dot(mmat.T,dat))-s
				deltas *= np.clip(np.max(np.abs(deltas)),None,0.5/steps[0])/np.max(np.abs(deltas))
				ds = 1-2*(c2p < chi2_th) # Direction sign; is chi squared too large or too small?
				c20 = np.mean((dat0-np.dot(Rij,np.exp(s+deltas*ds*steps[0])))**2.0/err**2.0)
				c21 = np.mean((dat0-np.dot(Rij,np.exp(s+deltas*ds*steps[1])))**2.0/err**2.0)
				interp_step = ((steps[0]*(c21-chi2_th)+steps[1]*(chi2_th-c20))/(c21-c20))
				s += deltas*ds*np.clip(interp_step,steps[0],steps[1])
				chi2[i,j] = np.mean((dat0-np.dot(Rij,np.exp(s)))**2/err**2)
				if((ds*(c2p-c20)/steps[0] < tol)*(k > kcon) or np.abs(chi2[i,j]-chi2_th) < tol): break
			dems[i,j,:] = np.exp(s)

	return dems,chi2