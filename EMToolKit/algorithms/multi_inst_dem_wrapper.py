"""
This module provides a wrapper for Differential Emission Measure (DEM) analysis
using multiple instruments. It includes functions for reprojecting data with
uncertainties and plotting the results.

Functions:
    multi_inst_simple_dem_wrapper: Wraps simple DEM analysis for multiple instruments.
    reproject_with_uncertainties: Reprojects a sequence of maps to the smallest map.
    plot_reprojected: Plots reprojected maps and their uncertainties.
"""

import matplotlib.pyplot as plt
import astropy.units as u
import os
import numpy as np

from ndcube import NDCube, NDCubeSequence, NDCollection
from sunpy.visualization.colormaps.color_tables import aia_color_table, xrt_color_table
from astropy.nddata import StdDevUncertainty
from EMToolKit.algorithms.simple_reg_dem_wrapper import autoloading_simple_reg_dem_wrapper

def multi_inst_simple_dem_wrapper(datasequence, wrapargs={}, doPlot=False, dat_dir="../data/multitest", recalc_simple=False):
    """
    Perform simple DEM analysis using multiple instruments.

    Parameters:
        datasequence (NDCubeSequence): Sequence of data cubes for multiple instruments.
        wrapargs (dict): Arguments for the DEM wrapper.
        doPlot (bool): Whether to plot the reprojected data immediately.
        dat_dir (str): Directory for data storage.
        recalc_simple (bool): Whether to recalculate simple DEM.

    Returns:
        NDCubeSequence: Resulting DEM analysis.
    """
    downprojected_sequence, nan_mask, coarsest_cube = reproject_with_uncertainties(datasequence)

    if doPlot:
        plot_reprojected(downprojected_sequence, nan_mask, coarsest_cube)

    # print("Performing simple DEM...")
    wrapargs["prepend"] = "multi_inst_"

    return autoloading_simple_reg_dem_wrapper(datasequence, dat_dir, wrapargs=wrapargs, recalc_simple=recalc_simple)


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
    for seq in datasequence:
        sz = seq.meta['CDELT1'] * seq.meta['CDELT2']
        if sz >= sz_min:
            sz_min = sz
            coarsest_cube = seq

    # Pull out the fine Sequence
    fine_sequence = NDCubeSequence([mp for mp in datasequence if mp is not coarsest_cube], meta=datasequence.meta)

    # Reproject the fine maps to the coarse map shape
    downprojected_sequence = NDCubeSequence([mp.reproject_to(coarsest_cube.wcs) for mp in fine_sequence])

    # Compute factor to scale the uncertainties by the area ratio
    orig_area = fine_sequence[0].data.shape[0] * fine_sequence[0].data.shape[1]
    new_area = downprojected_sequence[0].data.shape[0] * downprojected_sequence[0].data.shape[1]
    area_ratio_rt = np.sqrt(new_area / orig_area)

    # Reproject and scale the uncertainties
    uncertainties = [StdDevUncertainty(
        NDCube(area_ratio_rt * mp.uncertainty.array.data, mp.wcs, meta=mp.meta).
        reproject_to(coarsest_cube.wcs).data
    ) for mp in fine_sequence]
    for i in range(len(downprojected_sequence)):
        downprojected_sequence[i].uncertainty = uncertainties[i]

    # Combine the reprojected fine maps with the coarse map
    nan_mask = np.where(np.isnan(uncertainties[0].array), np.nan, 1)
    for cube in downprojected_sequence:
        nan_mask *= np.where(cube.data < nan_level, np.nan, 1)
    full_list = [x * nan_mask for x in downprojected_sequence]
    full_list.append(coarsest_cube * nan_mask)
    datasequence = NDCubeSequence(full_list, common_axis=0)

    print("Reprojected successfully!")

    return datasequence, nan_mask, coarsest_cube


def plot_reprojected(downprojected_sequence, nan_mask, coarsest_cube):
    """
    Plot the reprojected maps and their uncertainties.

    Parameters:
        downprojected_sequence (NDCubeSequence): Sequence of reprojected data cubes.
        nan_mask (numpy.ndarray): Mask for NaN values.
        coarsest_cube (NDCube): Coarsest data cube for comparison.
    """
    for aia_reproj_map in downprojected_sequence:
        fig = plt.figure(figsize=(12, 6))
        coarsest_cube *= nan_mask
        aia_reproj_map *= nan_mask
        ax1 = fig.add_subplot(1, 2, 1, projection=coarsest_cube)
        ax2 = fig.add_subplot(1, 2, 2, projection=coarsest_cube, sharex=ax1, sharey=ax1)
        ax1.set_facecolor("grey")
        ax2.set_facecolor("grey")
        ax1.imshow(coarsest_cube.data, cmap=xrt_color_table(), origin="lower")
        wave = aia_reproj_map.meta.get("wavelnth", 0) * u.angstrom
        # aia_reproj_map.plot( axes=ax1, alpha=0.75, cmap=aia_color_table(wave))
        ax1.imshow(np.sqrt(aia_reproj_map.data), alpha=0.75, cmap=aia_color_table(wave))
        ax2.imshow(aia_reproj_map.uncertainty.array, cmap="plasma")
        ax1.set_title(f'AIA {wave} overlaid on XRT\nshape: {aia_reproj_map.data.shape}')
        ax2.set_title(f'Uncertainty\nshape: {aia_reproj_map.uncertainty.array.shape}')
        plt.tight_layout()
        plt.show()