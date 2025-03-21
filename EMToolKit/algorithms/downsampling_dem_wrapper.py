"""
This is the multi-instrument capable Differential Emission Measure
(DEM) module using the *downsampling* method. It can handle input images
with differing spatial resolution (e.g., from different instruments such
as SDO/AIA and Hinode), and does so by by downsampling to the lowest common
denominator pixel grid of the inputs, then applying standard per-pixel
algorithms. It includes functions for reprojecting data with uncertainties
and plotting the results.

Functions:
    downsampling_dem_wrapper: Wraps simple DEM analysis for multiple instruments.
    reproject_with_uncertainties: Reprojects a sequence of maps to the smallest map.
    plot_reprojected: Plots reprojected maps and their uncertainties.
"""

import os, copy, numpy as np, astropy.units as u, matplotlib.pyplot as plt
from ndcube import NDCube, NDCubeSequence
from astropy.nddata import StdDevUncertainty

from ndcube import NDCube, NDCubeSequence, NDCollection
from sunpy.visualization.colormaps.color_tables import aia_color_table, xrt_color_table
from astropy.nddata import StdDevUncertainty
from EMToolKit.algorithms.simple_reg_dem_wrapper import autoloading_simple_reg_dem_wrapper
from EMToolKit.algorithms.sparse_em_wrapper import autoloading_sparse_em_wrapper


def downsampling_dem_wrapper(datasequence,*, wrapargs={}, method='simple', doPlot=False, save_dir=None, recalc=False):
    """
    Perform DEM analysis using multiple instruments. This is the main call for the downsampling method.

    Parameters:
        datasequence (NDCubeSequence): Sequence of data cubes for multiple instruments.
        wrapargs (dict): Arguments for the DEM wrapper.
        method (string): which DEM algorithm to use.
        doPlot (bool): Whether to plot the reprojected data immediately.
        dat_dir (str): Directory for data storage.
        recalc_simple (bool): Whether to recalculate simple DEM.

    Returns:
        NDCubeSequence: Resulting DEM analysis.
    """

    datasequence.meta = datasequence[0].meta

    downprojected_sequence, nan_mask, coarsest_cube = reproject_with_uncertainties(datasequence)

    if doPlot:
        plot_reprojected(downprojected_sequence, nan_mask, coarsest_cube)

    wrapargs["prepend"] = "multi_down_"
    if "simple" in method.casefold():
        return autoloading_simple_reg_dem_wrapper(downprojected_sequence, save_dir=save_dir, wrapargs=wrapargs, recalc=recalc)
    elif "sparse" in method.casefold():
        return autoloading_sparse_em_wrapper(downprojected_sequence, save_dir=save_dir, wrapargs=wrapargs, recalc=recalc)

def reproject_with_uncertainties(datasequence, nan_level=-50):
    """
    Reproject a sequence of maps to the smallest map with uncertainties.

    Parameters:
        datasequence (NDCubeSequence): Sequence of data cubes for multiple instruments.
        nan_level (float): Threshold level for NaN masking.

    Returns:
        tuple: Reprojected data sequence, NaN mask, and coarsest cube.
    """
    print("Reprojecting all maps to smallest map...", end="")
    # Find the map with the largest pixel scale
    sz_min = 0
    coarsest_cube = None
    for seq in datasequence:
        sz = seq.meta['CDELT1'] * seq.meta['CDELT2']
        if sz >= sz_min:
            sz_min = sz
            coarsest_cube = seq

    if coarsest_cube is None:
        raise ValueError("No coarsest cube found in the sequence")

    ## Pull out the fine Sequence
    #fine_sequence = NDCubeSequence([mp for mp in datasequence if mp is not coarsest_cube], meta=datasequence.meta)

    # Reproject the fine maps to the coarse map shape -- Note this DOES NOT pixel bin for large factors
    downprojected_sequence, uncertainties = [], []
    for mp in datasequence:
        mp_new = mp.reproject_to(coarsest_cube.wcs)
        err_copy = copy.deepcopy(mp)
        err_copy.data[:] = err_copy.uncertainty.array
        mp_new.uncertainty = StdDevUncertainty(err_copy.reproject_to(coarsest_cube.wcs).data)
        downprojected_sequence.append(mp_new)
        uncertainties.append(mp_new.uncertainty)
	#downprojected_sequence = NDCubeSequence([mp.reproject_to(coarsest_cube.wcs) for mp in datasequence], meta=datasequence.meta)
    #print(dir(downprojected_sequence[0]))

    ## Compute factor to scale the uncertainties by the area ratio (not apropos given the above)
    #orig_area = fine_sequence[0].data.shape[0] * fine_sequence[0].data.shape[1]
    #new_area = downprojected_sequence[0].data.shape[0] * downprojected_sequence[0].data.shape[1]
    area_ratio_rt = 1 #np.sqrt(new_area / orig_area)

    # print(f"Original area: {orig_area}, New area: {new_area}, Area ratio root: {area_ratio_rt}")

    # Reproject and scale the uncertainties
    #uncertainties = [StdDevUncertainty(
    #    NDCube(area_ratio_rt * mp.uncertainty.array.data, mp.wcs, meta=mp.meta).
    #    reproject_to(coarsest_cube.wcs).data
    #) for mp in downprojected_sequence]
    #for I in range(len(downprojected_sequence)):
    #    downprojected_sequence[I].uncertainty = uncertainties[I]

    ## Combine the reprojected fine maps with the coarse map
    nan_mask = np.isnan(uncertainties[0].array)
    for cube in downprojected_sequence:
        nan_mask *= np.isnan(cube.data) | (cube.data < nan_level)
    #full_list = [cube.data * nan_mask for cube in downprojected_sequence]
    #full_list.append(coarsest_cube.data * nan_mask)
    #datasequence = NDCubeSequence(full_list, common_axis=0, meta=datasequence.meta)

    print("Reprojected successfully!")

    return downprojected_sequence, nan_mask, coarsest_cube


import matplotlib.pyplot as plt
import numpy as np
import astropy.units as u
from ndcube import NDCubeSequence
from sunpy.map import Map

def plot_reprojected(downprojected_sequence, nan_mask, coarsest_cube):
    """
    Plot the reprojected maps and their uncertainties.

    Parameters:
        downprojected_sequence (NDCubeSequence): Sequence of reprojected data cubes.
        nan_mask (numpy.ndarray): Mask for NaN values.
        coarsest_cube (NDCube): Coarsest data cube for comparison.
    """
    for aia_reproj_map in Map(downprojected_sequence):
        fig = plt.figure(figsize=(12, 6))

        # Apply NaN mask to the data
        coarsest_data_masked = coarsest_cube.data * nan_mask
        aia_data_masked = aia_reproj_map.data * nan_mask
        # aia_uncertainty_masked = aia_reproj_map.uncertainty.array * nan_mask

        # Create axes with WCS projection
        ax1 = fig.add_subplot(1, 2, 1, projection=coarsest_cube.wcs)
        ax2 = fig.add_subplot(1, 2, 2, projection=coarsest_cube.wcs, sharex=ax1, sharey=ax1)

        ax1.set_facecolor("grey")
        ax2.set_facecolor("grey")

        # Plot the coarsest cube data
        ax1.imshow(coarsest_data_masked, cmap='gray', origin="lower")

        # Get wavelength for the AIA map
        wave = aia_reproj_map.meta.get("wavelnth", 0) * u.angstrom

        # Overlay the AIA reprojected data
        ax1.imshow(np.sqrt(aia_data_masked), alpha=0.75, cmap='inferno')

        # Plot the uncertainties
        # ax2.imshow(aia_uncertainty_masked, cmap="plasma")

        ax1.set_title(f'AIA {wave} overlaid on Coarsest Cube\nShape: {aia_reproj_map.data.shape}')
        ax2.set_title(f'Uncertainty\nShape: {aia_reproj_map.uncertainty.array.shape}')

        plt.tight_layout()
        plt.show()

# Helper function for color table (customize or replace with your actual function)
def aia_color_table(wave):
    return 'inferno'

# Example function call
# plot_reprojected(downprojected_sequence, nan_mask, coarsest_cube)
