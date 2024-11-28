import numpy as np
from sunpy.map import Map
from ndcube import NDCube, NDCubeSequence, NDCollection
from astropy.coordinates import SkyCoord
from astropy.nddata import StdDevUncertainty
from EMToolKit.algorithms.simple_reg_dem import simple_reg_dem
import os.path
import pickle
import time
import EMToolKit.EMToolKit as emtk


def simple_reg_dem_wrapper(datasequence, wrapargs=None):
    """
    Wrapper function for the simple regularized Differential Emission Measure (DEM) calculation.

    This function prepares input data and passes it to the `simple_reg_dem` algorithm.
    It processes the input data to ensure that all input maps have consistent dimensions,
    extracts the necessary metadata, and then calls the DEM calculation.

    Parameters
    ----------
    datasequence : NDCubeSequence
        A sequence of data cubes containing the observations. Each cube should contain
        2D spatial data with associated uncertainties and metadata.
    wrapargs : dict, optional
        Additional arguments to pass to the initialization routines of the `simple_reg_dem` function.

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
def simple_reg_dem_wrapper(datasequence, wrapargs=None):
    """
    Wrapper function for the simple regularized Differential Emission Measure (DEM) calculation.

    This function prepares input data and passes it to the `simple_reg_dem` algorithm.
    It processes the input data to ensure that all input maps have consistent dimensions,
    extracts the necessary metadata, and then calls the DEM calculation.

    Parameters
    ----------
    datasequence : NDCubeSequence
        A sequence of data cubes containing the observations. Each cube should contain
        2D spatial data with associated uncertainties and metadata.
    wrapargs : dict, optional
        Additional arguments to pass to the initialization routines of the `simple_reg_dem` function.

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

    logt = datasequence[0].meta['logt']  # Temperature bins
    datacube = np.zeros([nx, ny, nc])  # Data cube for storing observations
    errscube = np.zeros([nx, ny, nc])  # Data cube for storing uncertainties
    tresps = np.zeros([logt.size, nc])  # Temperature response functions
    exptimes = np.zeros(nc)  # Exposure times

    # Fill the data, uncertainty, and metadata arrays
    for i in range(nc):
        datacube[:, :, i] = datasequence[i].data[0:nx, 0:ny]
        errscube[:, :, i] = datasequence[i].uncertainty.array[0:nx, 0:ny]
        tresps[:, i] = datasequence[i].meta['tresp']
    [nx, ny] = datasequence[0].data.shape
    for seq in datasequence:
        [nx, ny] = [np.min([seq.data.shape[0], nx]), np.min([seq.data.shape[1], ny])]

    logt = datasequence[0].meta['logt']  # Temperature bins
    datacube = np.zeros([nx, ny, nc])  # Data cube for storing observations
    errscube = np.zeros([nx, ny, nc])  # Data cube for storing uncertainties
    tresps = np.zeros([logt.size, nc])  # Temperature response functions
    exptimes = np.zeros(nc)  # Exposure times

    # Fill the data, uncertainty, and metadata arrays
    for i in range(nc):
        datacube[:, :, i] = datasequence[i].data[0:nx, 0:ny]
        errscube[:, :, i] = datasequence[i].uncertainty.array[0:nx, 0:ny]
        tresps[:, i] = datasequence[i].meta['tresp']
        exptimes[i] = datasequence[i].meta['exptime']

    # Perform the DEM calculation
    coeffs, chi2 = simple_reg_dem(datacube, errscube, exptimes, logt, tresps)

    # Perform the DEM calculation
    coeffs, chi2 = simple_reg_dem(datacube, errscube, exptimes, logt, tresps)

    # Simple_reg_dem puts the temperature axis last. Transpose so it's the first:
    coeffs = coeffs.transpose(np.roll(np.arange(coeffs.ndim), 1))
    coeffs = coeffs.transpose(np.roll(np.arange(coeffs.ndim), 1))

    nt = logt.size
    wcs = datasequence[0].wcs  # WCS information from the first cube
    basislogt = np.linspace(np.min(logt), np.max(logt), 2 * (nt - 1) + 1)
    logts, bases = [], []

    # Create the basis functions for the temperature bins
    for i in range(nt):
    wcs = datasequence[0].wcs  # WCS information from the first cube
    basislogt = np.linspace(np.min(logt), np.max(logt), 2 * (nt - 1) + 1)
    logts, bases = [], []

    # Create the basis functions for the temperature bins
    for i in range(nt):
        basis = (basislogt == logt[i]).astype(np.float64)
        if i > 0:
            basis += (basislogt - logt[i - 1]) * (basislogt < logt[i]) * (basislogt > logt[i - 1]) / (logt[i] - logt[i - 1])
        if i < nt - 1:
            basis += (logt[i + 1] - basislogt) * (basislogt < logt[i + 1]) * (basislogt > logt[i]) / (logt[i + 1] - logt[i])
        if i > 0:
            basis += (basislogt - logt[i - 1]) * (basislogt < logt[i]) * (basislogt > logt[i - 1]) / (logt[i] - logt[i - 1])
        if i < nt - 1:
            basis += (logt[i + 1] - basislogt) * (basislogt < logt[i + 1]) * (basislogt > logt[i]) / (logt[i + 1] - logt[i])
        bases.append(basis)
        logts.append(basislogt)

    return list(coeffs), logts, bases, wcs, 'simple_reg_dem'


def autoloading_simple_reg_dem_wrapper(datasequence, data_dir=".data/default/", recalc=False, wrapargs={}):
    """
    Wrapper function that calculates or loads a precomputed simple regularized DEM.

    This function first checks if a precomputed DEM exists in the specified directory. If not,
    it calculates the DEM using `simple_reg_dem_wrapper`, saves the result, and returns it.

    Parameters
    ----------
    datasequence : NDCubeSequence
        A sequence of data cubes containing the observations. Each cube should contain
        2D spatial data with associated uncertainties and metadata.
    data_dir : str, optional
        The directory where the DEM result will be saved or loaded from. Default is ".data/default/".
    recalc_simple : bool, optional
        If True, the DEM will be recalculated even if a precomputed result exists. Default is False.
    wrapargs : dict, optional
        Additional arguments to pass to the initialization routines of the `simple_reg_dem_wrapper` function.

    Returns
    -------
    emtk.dem_model
        The DEM model generated from the input data.
    tuple
        The output from the `simple_reg_dem_wrapper` function.
    """
    # Create the directory if it does not exist
    if not os.path.exists(data_dir):
        os.makedirs(data_dir)

    pk_file = os.path.join(data_dir, wrapargs.get("prepend",'single_') + 'simple_reg_demsequence.pkl')

    # Load or calculate the DEM sequence
    if os.path.exists(pk_file) and not recalc:
        print('Loading simple_reg_demsequence from', pk_file)
        with open(pk_file, 'rb') as file:
            (simple_reg_demsequence, simpl_out) = pickle.load(file)
    else:
        print(f"Calculating {pk_file} from scratch...", end="")
        tstart = time.time()
        simpl_out = simple_reg_dem_wrapper(datasequence, wrapargs)
        print('Done! Simple method took', time.time() - tstart)
        print('Done! Simple method took', time.time() - tstart)
        simple_reg_demsequence = emtk.dem_model(*simpl_out, simple_reg_dem_wrapper)

        # Save the DEM sequence to a file

        # Save the DEM sequence to a file
        with open(pk_file, 'wb') as file:
            pickle.dump((simple_reg_demsequence, simpl_out), file)


    return simple_reg_demsequence, simpl_out