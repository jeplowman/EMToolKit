"""
This module provides functions to initialize and solve simple regularized differential emission measure (DEM)
inversions using a sequence of solar images. The primary functions are `simple_reg_dem_wrapper` and
`autoloading_simple_reg_dem_wrapper`.

Functions:
    - simple_reg_dem_wrapper: Wraps the simple_reg_dem initialization and solving routines.
    - autoloading_simple_reg_dem_wrapper: Loads or recalculates the simple_reg_dem DEM sequence and saves the result to a pickle file.

Usage:
    This module can be used to process a sequence of solar images and compute the regularized DEM
    using the `simple_reg_dem_wrapper` function. The `autoloading_simple_reg_dem_wrapper` function provides
    an option to cache the results for faster loading.

Example:
    import numpy as np
    from sunpy.map import Map
    from ndcube import NDCube, NDCubeSequence, NDCollection
    from astropy.coordinates import SkyCoord
    from astropy.nddata import StdDevUncertainty
    from EMToolKit.algorithms.simple_reg_dem import simple_reg_dem
    import os
    import pickle
    import time
    import EMToolKit.EMToolKit as emtk

    datasequence = ...  # Load or create your data sequence
    simple_reg_demsequence, simpl_out = autoloading_simple_reg_dem_wrapper(datasequence)
"""

import numpy as np
from sunpy.map import Map
from ndcube import NDCube, NDCubeSequence, NDCollection
from astropy.coordinates import SkyCoord
from astropy.nddata import StdDevUncertainty
from EMToolKit.algorithms.simple_reg_dem import simple_reg_dem
import os
import pickle
import time
import EMToolKit.EMToolKit as emtk

def simple_reg_dem_wrapper(datasequence, wrapargs=None):
    """
    Wraps the simple_reg_dem initialization and solving routines.

    Parameters:
    datasequence (NDCubeSequence): Input data sequence.
    wrapargs (dict): Additional arguments for simple_reg_dem routines.

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
    tresps = np.zeros((logt.size, nc))
    exptimes = np.zeros(nc)

    for i in range(nc):
        datacube[:, :, i] = datasequence[i].data[:nx, :ny]
        # errscube[:, :, i] = datasequence[i].uncertainty.array[:nx, :ny]
        tresps[:, i] = datasequence[i].meta['tresp']
        exptimes[i] = datasequence[i].meta['exptime']

    coeffs, chi2 = simple_reg_dem(datacube, datacube, exptimes, logt, tresps)
    # Simple_reg_dem puts the temperature axis last. Transpose so it's the first:
    coeffs = coeffs.transpose(np.roll(np.arange(coeffs.ndim), 1))

    nt = logt.size
    wcs = datasequence[0].wcs
    basislogt = np.linspace(np.min(logt), np.max(logt), 2 * (nt - 1) + 1)
    logts, bases = [], []

    for i in range(nt):
        basis = (basislogt == logt[i]).astype(np.float64)
        if i > 0:
            basis += (basislogt - logt[i - 1]) * (basislogt < logt[i]) * (basislogt > logt[i - 1]) / (logt[i] - logt[i - 1])
        if i < nt - 1:
            basis += (logt[i + 1] - basislogt) * (basislogt < logt[i + 1]) * (basislogt > logt[i]) / (logt[i + 1] - logt[i])
        bases.append(basis)
        logts.append(basislogt)

    return list(coeffs), logts, bases, wcs, 'simple_reg_dem'

def autoloading_simple_reg_dem_wrapper(datasequence, data_dir="../", recalc_simple=False, wrapargs=None):
    """
    Loads or recalculates the simple_reg_dem DEM sequence and saves the result to a pickle file.

    Parameters:
    datasequence (NDCubeSequence): Input data sequence.
    data_dir (str): Directory to save or load the pickle file.
    recalc_simple (bool): Force recalculation of simple_reg_dem if True.
    wrapargs (dict): Additional arguments for simple_reg_dem routines.

    Returns:
    tuple: simple_reg_demsequence and simpl_out.
    """
    pk_file = os.path.join(data_dir, 'simple_reg_demsequence.pkl')

    if os.path.exists(pk_file) and not recalc_simple:
        print('Loading simple_reg_demsequence from', pk_file)
        try:
            with open(pk_file, 'rb') as file:
                simple_reg_demsequence, simpl_out = pickle.load(file)
        except (OSError, IOError, pickle.UnpicklingError) as e:
            print(f"Error loading pickle file: {e}")
            recalc_simple = True
    else:
        print("Calculating DEM from scratch...", end="")
        tstart = time.time()
        simpl_out = simple_reg_dem_wrapper(datasequence, wrapargs)
        print('Done! Simple method took', time.time() - tstart)
        simple_reg_demsequence = emtk.dem_model(*simpl_out, simple_reg_dem_wrapper)
        if not os.path.exists(data_dir):
            os.makedirs(data_dir)
        try:
            with open(pk_file, 'wb') as file:
                pickle.dump((simple_reg_demsequence, simpl_out), file)
        except (OSError, IOError, pickle.PicklingError) as e:
            print(f"Error saving pickle file: {e}")

    return simple_reg_demsequence, simpl_out