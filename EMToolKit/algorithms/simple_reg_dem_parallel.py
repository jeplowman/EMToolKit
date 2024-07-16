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

import numpy as np
from scipy.linalg import cho_factor, cho_solve
from concurrent.futures import ProcessPoolExecutor, as_completed
from tqdm import tqdm

def precompute_matrices(data, errors, exptimes, logt, tresps, drv_con):
    """
    Precomputes matrices and vectors required for regularized differential emission measure.

    Parameters:
    - data: The observed data array.
    - errors: The uncertainties associated with the data.
    - exptimes: Exposure times for each observation.
    - logt: Logarithm of temperature bins.
    - tresps: Temperature response functions.
    - drv_con: Drive constant for regularization.

    Returns:
    - nt_ones: Array of ones of size nt.
    - nx, ny: Dimensions of the data array.
    - Bij: Matrix used in regularization.
    - Rij: Matrix mapping coefficients to data.
    - Dij: Matrix used in regularization.
    - regmat: Regularization matrix.
    - rvec: Sum of Rij across axis=1.
    """
    nt, nd = tresps.shape
    nt_ones = np.ones(nt)
    nx, ny, _ = data.shape
    dT = logt[1:nt] - logt[0:nt-1]
    dTleft, dTright = np.diag(np.hstack([dT, 0])), np.diag(np.hstack([0, dT]))
    idTleft, idTright = np.diag(np.hstack([1.0/dT, 0])), np.diag(np.hstack([0, 1.0/dT]))
    Bij = ((dTleft + dTright) * 2.0 + np.roll(dTright, -1, axis=0) + np.roll(dTleft, 1, axis=0)) / 6.0
    Rij = np.matmul((tresps * np.outer(nt_ones, exptimes)).T, Bij)  # Matrix mapping coefficients to data
    Dij = idTleft + idTright - np.roll(idTright, -1, axis=0) - np.roll(idTleft, 1, axis=0)
    regmat = Dij * nd / (drv_con ** 2 * (logt[nt - 1] - logt[0]))
    rvec = np.sum(Rij, axis=1)

    return nt_ones, nx, ny, nt, Bij, Rij, Dij, regmat, rvec



def process_pixel(i, j, data, errors, Rij, regmat, rvec, nt_ones, kmax, kcon, steps, drv_con, chi2_th, tol):
    """
    Processes a single pixel (i, j) in a dataset to compute its differential emission measure (DEM) and chi-squared value.

    Parameters:
    - i, j: Indices of the pixel in the data array.
    - data: The observed data array.
    - errors: The uncertainties associated with the data.
    - Rij: Matrix mapping coefficients to data.
    - regmat: Regularization matrix.
    - rvec: Sum of Rij across axis=1.
    - nt_ones: Array of ones of size nt.
    - kmax: Maximum number of iterations.
    - kcon: Convergence criterion based on the number of iterations.
    - steps: Step sizes for the iterative method.
    - drv_con: Drive constant for regularization.
    - chi2_th: Chi-squared threshold for stopping criterion.
    - tol: Tolerance for convergence criterion.

    Returns:
    - i, j: Indices of the processed pixel.
    - dem_pixel: Computed DEM value for the pixel.
    - chi2_pixel: Computed chi-squared value for the pixel.

    Notes:
    - The function skips pixels with NaN values in data or errors.
    - If the Cholesky factorization fails, the iteration breaks early.
    - The function employs an iterative method to adjust the DEM based on the observed data, errors, and a regularization term.
    """
    dat00 = data[i, j, :]
    nan_level = -10 # This is the level below which we will mask out the data
    if np.isnan(dat00).any() or np.any(dat00 < nan_level):
        # print(f'Pix: Skipping pixel {i}, {j} with NaNs in data or errors')
        return i, j, np.nan, np.nan

    err = errors[i, j, :]
    dat0 = np.clip(data[i, j, :], 0.0, None)
    s = np.log(np.sum((rvec) * (np.clip(dat0,1.0e-2,None) / err**2)) / np.sum((rvec / err)**2) / nt_ones)
    for k in range(kmax):
        dat = (dat0 - np.matmul(Rij, ((1 - s) * np.exp(s)))) / err
        mmat = Rij * np.outer(1.0 / err, np.exp(s))
        amat = np.matmul(mmat.T, mmat) + regmat
        try:
            c, low = cho_factor(amat)
        except np.linalg.LinAlgError:
            break
        c2p = np.mean((dat0 - np.dot(Rij, np.exp(s)))**2 / err**2)
        deltas = cho_solve((c, low), np.dot(mmat.T, dat)) - s
        deltas *= np.clip(np.max(np.abs(deltas)), None, 0.5 / steps[0]) / np.max(np.abs(deltas))
        ds = 1 - 2 * (c2p < chi2_th)  # Direction sign
        c20 = np.mean((dat0 - np.dot(Rij, np.exp(s + deltas * ds * steps[0])))**2.0 / err**2.0)
        c21 = np.mean((dat0 - np.dot(Rij, np.exp(s + deltas * ds * steps[1])))**2.0 / err**2.0)
        interp_step = ((steps[0] * (c21 - chi2_th) + steps[1] * (chi2_th - c20)) / (c21 - c20))
        s += deltas * ds * np.clip(interp_step, steps[0], steps[1])
        chi2_pixel = np.mean((dat0 - np.dot(Rij, np.exp(s)))**2 / err**2)
        if (ds * (c2p - c20) / steps[0] < tol) * (k > kcon) or np.abs(chi2_pixel - chi2_th) < tol:
            break
    dem_pixel = np.exp(s)
    return i, j, dem_pixel, chi2_pixel



def simple_reg_dem_serial(data, errors, exptimes, logt, tresps,
                          kmax=100, kcon=5, steps=[0.1,0.5], drv_con=8.0, chi2_th=1.0, tol=0.1):
    """
    Perform a serial computation of differential emission measure (DEM) and chi-squared values
    for observational data against a set of response functions.

    Parameters:
    - data: Observed data, a 3D array (nx, ny, nd).
    - errors: Uncertainties in the observed data, same shape as data.
    - exptimes: Exposure times for each observation, 1D array of length nd.
    - logt: Logarithm of temperature bins, 1D array.
    - tresps: Temperature response functions, 2D array (nt, nd).
    - kmax: Maximum number of iterations for the fitting process.
    - kcon: Minimum number of iterations before convergence can be considered.
    - steps: Step sizes to use in the fitting iterations, list of two floats.
    - drv_con: Drive constant for the regularization term.
    - chi2_th: Threshold chi-squared value for stopping the iterations.
    - tol: Tolerance for convergence.

    Returns:
    - dems: Computed DEM values, a 3D array (nx, ny, nt).
    - chi2: Computed chi-squared values for each pixel, a 2D array (nx, ny).
    """
    # Precompute matrices based on the input parameters
    nt_ones, nx, ny, nt, Bij, Rij, Dij, regmat, rvec = precompute_matrices(data, errors, exptimes, logt, tresps, drv_con)

    # Initialize the arrays to store DEMs and chi-squared values
    dems = np.zeros([nx,ny,nt])
    chi2 = np.zeros([nx,ny]) - 1.0

    # Create a progress bar to monitor the computation
    total_pixels = nx * ny
    with tqdm(total=total_pixels, desc="\tProcessing Pixels (Serial)") as pbar:
        # Iterate over each pixel to compute its DEM and chi-squared value
        for i in range(nx):
            for j in range(ny):
                i, j, dems_pixel, chi2_pixel = process_pixel(i, j, data, errors, Rij, regmat, rvec,
                                                             nt_ones, kmax, kcon, steps, drv_con, chi2_th, tol)
                dems[i, j, :] = dems_pixel
                chi2[i, j] = chi2_pixel
                pbar.update(1)  # Update the progress bar for each processed pixel

    return dems, chi2



def simple_reg_dem_parallel(data, errors, exptimes, logt, tresps,
                            kmax=100, kcon=5, steps=[0.1,0.5], drv_con=8.0, chi2_th=1.0, tol=0.1):
    """
    Perform a parallel computation of differential emission measure (DEM) and chi-squared values
    for observational data against a set of response functions using multiple processes.

    Parameters are the same as in simple_reg_dem_serial.

    Returns:
    - dems: Computed DEM values, a 3D array (nx, ny, nt).
    - chi2: Computed chi-squared values for each pixel, a 2D array (nx, ny).
    """
    # Precompute matrices based on the input parameters
    nt_ones, nx, ny, nt, Bij, Rij, Dij, regmat, rvec = precompute_matrices(data, errors, exptimes, logt, tresps, drv_con)

    # Initialize the arrays to store DEMs and chi-squared values
    dems = np.zeros([nx, ny, nt])
    chi2 = np.zeros([nx, ny]) - 1.0

    # Use ProcessPoolExecutor for executing computations in parallel
    with ProcessPoolExecutor() as executor:
        # Create all futures for processing each pixel in parallel
        futures = {executor.submit(process_pixel, i, j, data, errors, Rij, regmat, rvec, nt_ones,
                                   kmax, kcon, steps, drv_con, chi2_th, tol): (i, j) for i in range(nx) for j in range(ny)}

        # Track progress of futures as they complete
        for future in tqdm(as_completed(futures), total=len(futures), desc="\tProcessing Pixels (Parallel)"):
            i, j, dems_pixel, chi2_pixel = future.result()
            dems[i, j, :] = dems_pixel
            chi2[i, j] = chi2_pixel

    return dems, chi2



def simple_reg_dem(*args, doParallel=True, **kwargs):
    """
    Wrapper function to call the serial or parallel version of the simple regularized DEM computation.

    Accepts the same parameters as simple_reg_dem_serial and simple_reg_dem_parallel.

    Returns:
    - dems: Computed DEM values, a 3D array (nx, ny, nt).
    - chi2: Computed chi-squared values for each pixel, a 2D array (nx, ny).
    """
    if doParallel:
        return simple_reg_dem_parallel(*args, **kwargs)
    else:
        return simple_reg_dem_serial(*args, **kwargs)
