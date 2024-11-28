import copy
import numpy as np
from sunpy.map import Map
from ndcube import NDCube, NDCubeSequence, NDCollection
from astropy.coordinates import SkyCoord
from astropy.nddata import StdDevUncertainty
from EMToolKit.algorithms.sparse_em import sparse_em_init, sparse_em_solve
import os.path
import pickle
import time
import EMToolKit.EMToolKit as emtk


def sparse_em_wrapper(datasequence, wrapargs={}):
    """
    Wrapper function for the sparse regularized Differential Emission Measure (DEM) calculation.

    This function prepares input data and passes it to the `sparse_em_init` and `sparse_em_solve`
    algorithms. It processes the input data to ensure that all input maps have consistent dimensions,
    extracts the necessary metadata, and then calls the DEM calculation.

    Parameters
    ----------
    datasequence : NDCubeSequence
        A sequence of data cubes containing the observations. Each cube should contain
        2D spatial data with associated uncertainties and metadata.
    wrapargs : dict, optional
        Additional arguments to pass to the initialization routines of the `sparse_em` functions.

    Returns
    -------
    list of numpy.ndarray
        The calculated DEM coefficients for each temperature bin.
    list of numpy.ndarray
        The temperature bins used in the calculation.
    list of numpy.ndarray
        The basis functions for the temperature bins.
    WCS
        World Coordinate System (WCS) information from the input data.
    str
        A string identifier for the DEM method used.
    """
    # Determine the smallest data cube dimensions in the sequence
    nc = len(datasequence)
    [nx, ny] = datasequence[0].data.shape
    for seq in datasequence:
        [nx, ny] = [np.min([seq.data.shape[0], nx]), np.min([seq.data.shape[1], ny])]

    # logt = datasequence[0].meta['logt']  # Temperature bins
    datacube = np.zeros([nx, ny, nc])  # Data cube for storing observations
    errscube = np.zeros([nx, ny, nc])  # Data cube for storing uncertainties
    trlogts, tresps = [], []  # Temperature response functions and temperature bins
    exptimes = np.zeros(nc)  # Exposure times

    # Fill the data, uncertainty, and metadata arrays
    for i in range(nc):
        datacube[:, :, i] = datasequence[i].data[0:nx, 0:ny]
        errscube[:, :, i] = datasequence[i].uncertainty.array[0:nx, 0:ny]
        tresps.append(datasequence[i].meta['tresp'])
        trlogts.append(datasequence[i].meta['logt'])
        exptimes[i] = datasequence[i].meta['exptime']

    # Initialize and solve the sparse EM problem
    Dict, lgtaxis, basis_funcs, bases_sigmas = sparse_em_init(trlogts, tresps, differential=True)
    coeffs, zmax, status = sparse_em_solve(datacube, errscube, exptimes, Dict)

    # Sparse_em_solve puts the temperature axis last. Transpose so it's the first:
    coeffs = coeffs.transpose(np.roll(np.arange(coeffs.ndim), 1))
    # Initialize and solve the sparse EM problem
    Dict, lgtaxis, basis_funcs, bases_sigmas = sparse_em_init(trlogts, tresps, differential=True)
    coeffs, zmax, status = sparse_em_solve(datacube, errscube, exptimes, Dict)

    # Sparse_em_solve puts the temperature axis last. Transpose so it's the first:
    coeffs = coeffs.transpose(np.roll(np.arange(coeffs.ndim), 1))

    nchannels, nb = Dict.shape
    wcs = datasequence[0].wcs  # WCS information from the first cube
    logts = nb * [lgtaxis]

    name = wrapargs.get("prepend", 'single_') + 'sparse_em'

    return list(coeffs), logts, list(basis_funcs), wcs, name



def autoloading_sparse_em_wrapper(datasequence, data_dir=".data/default", *, wrapargs={}, recalc=False):
    """
    Wrapper function that calculates or loads a precomputed sparse regularized DEM.

    This function first checks if a precomputed DEM exists in the specified directory. If not,
    it calculates the DEM using `sparse_em_wrapper`, saves the result, and returns it.

    Parameters
    ----------
    datasequence : NDCubeSequence
        A sequence of data cubes containing the observations. Each cube should contain
        2D spatial data with associated uncertainties and metadata.
    data_dir : str, optional
        The directory where the DEM result will be saved or loaded from. Default is "../.data/default".
    recalc : bool, optional
        If True, the DEM will be recalculated even if a precomputed result exists. Default is False.
    wrapargs : A dictionary

    Returns
    -------
    emtk.dem_model
        The DEM model generated from the input data.
    tuple
        The output from the `sparse_em_wrapper` function.
    """
    # Create the directory if it does not exist
    if not os.path.exists(data_dir):
        os.makedirs(data_dir)

    pk_file = os.path.join(data_dir, wrapargs.get("prepend",'single_') + 'sparse_em_demsequence.pkl')

    # Load or calculate the DEM sequence
    if os.path.exists(pk_file) and not recalc:
        print(f'Loading sparse_em_demsequence from {pk_file}')
        with open(pk_file, 'rb') as file:
            sparse_em_demsequence, spars_out = pickle.load(file)
    else:
        tstart = time.time()
        print(f"Calculating {pk_file} from scratch...", end="")

        spars_out = sparse_em_wrapper(datasequence, wrapargs=wrapargs)
        sparse_em_demsequence = emtk.dem_model(*spars_out, sparse_em_wrapper)
        print('Sparse method took', time.time() - tstart)

        # Save the DEM sequence to a file
        with open(pk_file, 'wb') as file:
            pickle.dump((sparse_em_demsequence, spars_out), file)

    return sparse_em_demsequence, spars_out

