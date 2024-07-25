"""
This module provides functions to initialize and solve sparse emission measure (EM) inversions using
a sequence of solar images. The primary functions are `sparse_em_wrapper` and `autoloading_sparse_em_wrapper`.

Functions:
    - sparse_em_wrapper: Wraps the sparse_em initialization and solving routines.
    - autoloading_sparse_em_wrapper: Loads or recalculates the sparse_em DEM sequence and saves the result to a pickle file.

Usage:
    This module can be used to process a sequence of solar images and compute the sparse emission measure
    using the `sparse_em_wrapper` function. The `autoloading_sparse_em_wrapper` function provides an option
    to cache the results for faster loading.

Example:
    import copy
    import numpy as np
    from sunpy.map import Map
    from ndcube import NDCube, NDCubeSequence, NDCollection
    from astropy.coordinates import SkyCoord
    from astropy.nddata import StdDevUncertainty
    from EMToolKit.algorithms.sparse_em import sparse_em_init, sparse_em_solve
    import os
    import pickle
    import time
    import EMToolKit.EMToolKit as emtk

    datasequence = ...  # Load or create your data sequence
    sparse_em_demsequence, spars_out = autoloading_sparse_em_wrapper(datasequence)
"""

import copy
import numpy as np
from sunpy.map import Map
from ndcube import NDCube, NDCubeSequence, NDCollection
from astropy.coordinates import SkyCoord
from astropy.nddata import StdDevUncertainty
from EMToolKit.algorithms.sparse_em import sparse_em_init, sparse_em_solve
import os
import pickle
import time
import EMToolKit.EMToolKit as emtk

def sparse_em_wrapper(datasequence, wrapargs=None):
    """
    Wraps the sparse_em initialization and solving routines.

    Parameters:
    datasequence (NDCubeSequence): Input data sequence.
    wrapargs (dict): Additional arguments for sparse_em routines.

    Returns:
    tuple: Coefficients, log temperatures, basis functions, WCS, and method name.
    """
    nc = len(datasequence)
    nx, ny = datasequence[0].data.shape
    for seq in datasequence:
        nx, ny = min(seq.data.shape[0], nx), min(seq.data.shape[1], ny)

    logt = datasequence[0].meta['logt']
    datacube = np.zeros((nx, ny, nc))
    errscube = np.zeros((nx, ny, nc))
    tresps = [seq.meta['tresp'] for seq in datasequence]
    trlogts = [seq.meta['logt'] for seq in datasequence]
    exptimes = np.array([seq.meta['exptime'] for seq in datasequence])

    for i in range(nc):
        datacube[:, :, i] = datasequence[i].data[:nx, :ny]
        errscube[:, :, i] = datasequence[i].uncertainty.array[:nx, :ny]

    Dict, lgtaxis, basis_funcs, bases_sigmas = sparse_em_init(trlogts, tresps, differential=True)
    coeffs, zmax, status = sparse_em_solve(datacube, datacube, exptimes, Dict)
    coeffs = coeffs.transpose(np.roll(np.arange(coeffs.ndim), 1))

    nchannels, nb = Dict.shape
    wcs = datasequence[0].wcs
    logts = [lgtaxis] * nb

    return list(coeffs), logts, list(basis_funcs), wcs, 'sparse_em'

def autoloading_sparse_em_wrapper(datasequence, data_dir="../data/multitest", recalc_sparse=False, wrapargs=None):
    """
    Loads or recalculates the sparse_em DEM sequence and saves the result to a pickle file.

    Parameters:
    datasequence (NDCubeSequence): Input data sequence.
    data_dir (str): Directory to save or load the pickle file.
    recalc_sparse (bool): Force recalculation of sparse_em if True.
    wrapargs (dict): Additional arguments for sparse_em routines.

    Returns:
    tuple: sparse_em_demsequence and spars_out.
    """
    pk_file = os.path.join(data_dir, 'sparse_em_demsequence.pkl')

    if os.path.exists(pk_file) and not recalc_sparse:
        try:
            with open(pk_file, 'rb') as file:
                sparse_em_demsequence, spars_out = pickle.load(file)
        except (OSError, IOError, pickle.UnpicklingError) as e:
            print(f"Error loading pickle file: {e}")
            recalc_sparse = True

    if not os.path.exists(pk_file) or recalc_sparse:
        tstart = time.time()
        spars_out = sparse_em_wrapper(datasequence, wrapargs=wrapargs)
        sparse_em_demsequence = emtk.dem_model(*spars_out, sparse_em_wrapper)
        print('Sparse method took', time.time() - tstart)
        try:
            with open(pk_file, 'wb') as file:
                pickle.dump((sparse_em_demsequence, spars_out), file)
        except (OSError, IOError, pickle.PicklingError) as e:
            print(f"Error saving pickle file: {e}")

    return sparse_em_demsequence, spars_out