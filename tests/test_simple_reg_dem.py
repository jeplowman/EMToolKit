"""
Tests that the Numba implementation of simple_reg_dem reproduces the original
NumPy/SciPy one, which is kept here as the reference.
"""

import numba
import numpy as np
import pytest
from scipy.linalg import cho_factor, cho_solve

from EMToolKit.algorithms.simple_reg_dem import simple_reg_dem
from EMToolKit.instruments.aia import AIA_TEMPERATURE_RESPONSE_TABLE, AIA_TEMPERATURES


def simple_reg_dem_original(data: np.ndarray, errors: np.ndarray, exptimes: np.ndarray, logt: np.ndarray, tresps: np.ndarray,
                            kmax: int = 100, kcon: int = 5, steps: list[float] = [0.1, 0.5], drv_con: float = 8.0,
                            chi2_th: float = 1.0, tol: float = 0.1) -> tuple[np.ndarray, np.ndarray]:
    """simple_reg_dem as it was before it was compiled with Numba, line for line."""
    [nt, nd] = tresps.shape
    nt_ones = np.ones(nt)
    [nx, ny, nd] = data.shape
    dT = logt[1:nt]-logt[0:nt-1]
    [dTleft, dTright] = [np.diag(np.hstack([dT, 0])), np.diag(np.hstack([0, dT]))]
    [idTleft, idTright] = [np.diag(np.hstack([1.0/dT, 0])), np.diag(np.hstack([0, 1.0/dT]))]
    Bij = ((dTleft+dTright)*2.0 + np.roll(dTright, -1, axis=0) + np.roll(dTleft, 1, axis=0))/6.0
    Rij = np.matmul((tresps*np.outer(nt_ones, exptimes)).T, Bij)  # Matrix mapping coefficents to data
    Dij = idTleft+idTright - np.roll(idTright, -1, axis=0) - np.roll(idTleft, 1, axis=0)
    regmat = Dij*nd/(drv_con**2*(logt[nt-1]-logt[0]))
    rvec = np.sum(Rij, axis=1)

    dems = np.zeros([nx, ny, nt])
    chi2 = np.zeros([nx, ny]) - 1.0
    for i in range(0, nx):
        for j in range(0, ny):
            err = errors[i, j, :]
            dat0 = np.clip(data[i, j, :], 0.0, None)
            s = np.log(np.sum((rvec)*((np.clip(dat0, 1.0e-2, None))/err**2))/np.sum((rvec/err)**2)/nt_ones)
            for k in range(0, kmax):
                dat = (dat0-np.matmul(Rij, ((1-s)*np.exp(s))))/err  # Correct data by f(s)-s*f'(s)...
                mmat = Rij*np.outer(1.0/err, np.exp(s))
                amat = np.matmul(mmat.T, mmat)+regmat
                try: [c, low] = cho_factor(amat)
                except: break
                c2p = np.mean((dat0-np.dot(Rij, np.exp(s)))**2/err**2)
                deltas = cho_solve((c, low), np.dot(mmat.T, dat))-s
                deltas *= np.clip(np.max(np.abs(deltas)), None, 0.5/steps[0])/np.max(np.abs(deltas))
                ds = 1-2*(c2p < chi2_th)  # Direction sign; is chi squared too large or too small?
                c20 = np.mean((dat0-np.dot(Rij, np.exp(s+deltas*ds*steps[0])))**2.0/err**2.0)
                c21 = np.mean((dat0-np.dot(Rij, np.exp(s+deltas*ds*steps[1])))**2.0/err**2.0)
                interp_step = ((steps[0]*(c21-chi2_th)+steps[1]*(chi2_th-c20))/(c21-c20))
                s += deltas*ds*np.clip(interp_step, steps[0], steps[1])
                chi2[i, j] = np.mean((dat0-np.dot(Rij, np.exp(s)))**2/err**2)
                if ((ds*(c2p-c20)/steps[0] < tol)*(k > kcon) or np.abs(chi2[i, j]-chi2_th) < tol): break
            dems[i, j, :] = np.exp(s)

    return dems, chi2


logt = np.arange(5.5, 7.5, 0.05)
tresps = np.stack([np.interp(logt, AIA_TEMPERATURES, r) for r in AIA_TEMPERATURE_RESPONSE_TABLE.T], axis=-1)
exptimes = np.array([2.9, 2.9, 2.0, 2.0, 2.9, 2.9])
gains = np.array([2.128, 1.523, 1.168, 1.024, 0.946, 0.596])  # DN per photon, as in estimate_aia_error


def observation(nx: int, ny: int, seed: int = 0) -> tuple[np.ndarray, np.ndarray]:
    """
    Counts and their uncertainties, (nx, ny, channel), from the random DEMs of Section 3.1
    of Plowman & Caspi (2020): sums of five Gaussians in log T, with total emission measures
    from 5e27 to 5e28 cm^-5, widths 0.05 to 0.15, and peaks from log T = 6 to 7, with
    photon noise and a read noise of 1.15 DN.
    """
    rng = np.random.default_rng(seed)
    em = 10**rng.uniform(np.log10(5e27), np.log10(5e28), (nx, ny, 5, 1))
    width = rng.uniform(0.05, 0.15, (nx, ny, 5, 1))
    center = rng.uniform(6.0, 7.0, (nx, ny, 5, 1))
    dem = (em/(np.sqrt(2*np.pi)*width)*np.exp(-0.5*((logt-center)/width)**2)).sum(-2)
    counts = (dem[..., None]*tresps).sum(-2)*(logt[1]-logt[0])*exptimes
    errors = np.sqrt(counts*gains + 1.15**2)
    return counts + rng.normal(size=counts.shape)*errors, errors


@pytest.mark.parametrize("kwargs", [
    dict(),
    dict(kcon=15, steps=[0.1, 0.75], drv_con=4.0),
    dict(chi2_th=2.0, tol=0.01),
])
def test_matches_original(kwargs: dict) -> None:
    """
    The DEMs and chi^2 agree with the original to rounding error. Not to the last bit,
    since the Cholesky factorizations round differently and the iteration amplifies it,
    but to about 1e-11.
    """
    data, errors = observation(8, 10)
    dems, chi2 = simple_reg_dem(data, errors, exptimes, logt, tresps, **kwargs)
    dems_original, chi2_original = simple_reg_dem_original(data, errors, exptimes, logt, tresps, **kwargs)
    assert dems.shape == dems_original.shape and chi2.shape == chi2_original.shape
    assert np.allclose(dems, dems_original, rtol=1e-6, atol=0)
    assert np.allclose(chi2, chi2_original, rtol=1e-6, atol=0)


def test_failed_pixels() -> None:
    """
    Pixels with a NaN datum, or a NaN or zero error, fail as in the original: a NaN DEM and
    chi^2 of -1. A negative error is taken as its magnitude, as in the original. The other
    pixels are not affected.
    """
    data, errors = observation(1, 6)
    data[0, 1, 2] = np.nan
    errors[0, 2, 3] = np.nan
    errors[0, 3, 4] = 0.0
    errors[0, 4, 5] *= -1
    dems, chi2 = simple_reg_dem(data, errors, exptimes, logt, tresps)
    with np.errstate(divide="ignore", invalid="ignore"):
        dems_original, chi2_original = simple_reg_dem_original(data, errors, exptimes, logt, tresps)
    assert np.array_equal(chi2[0, 1:4], [-1, -1, -1])
    assert np.array_equal(chi2 == -1, chi2_original == -1)
    assert np.array_equal(np.isnan(dems), np.isnan(dems_original))
    assert np.allclose(dems, dems_original, rtol=1e-6, atol=0, equal_nan=True)
    assert np.allclose(chi2, chi2_original, rtol=1e-6, atol=0)


def test_threads() -> None:
    """Each pixel gets the same result whichever thread runs it."""
    data, errors = observation(7, 9)
    threads = numba.get_num_threads()
    try:
        numba.set_num_threads(1)
        dems_one, chi2_one = simple_reg_dem(data, errors, exptimes, logt, tresps)
    finally:
        numba.set_num_threads(threads)
    dems, chi2 = simple_reg_dem(data, errors, exptimes, logt, tresps)
    assert np.array_equal(dems, dems_one) and np.array_equal(chi2, chi2_one)
