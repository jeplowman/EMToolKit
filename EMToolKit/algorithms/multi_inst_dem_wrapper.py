"""
This module provides a wrapper for Differential Emission Measure (DEM) analysis
using multiple instruments. It includes functions for reprojecting data with
uncertainties and plotting the results.

Functions:
    - multi_inst_simple_dem_wrapper: Wraps simple DEM analysis for multiple instruments.
    - find_coarsest_cube: Finds the cube with the largest pixel scale in the sequence.
    - reproject_cubes: Reprojects fine maps to the coarse map shape.
    - scale_uncertainties: Scales and reprojects uncertainties.
    - combine_sequences: Combines reprojected fine maps with the coarse map and applies NaN mask.
    - reproject_with_uncertainties: Reprojects a sequence of maps to the smallest map with uncertainties.
    - plot_reprojected: Plots reprojected maps and their uncertainties.
"""

import os
import numpy as np
import matplotlib.pyplot as plt
import astropy.units as u
from ndcube import NDCube, NDCubeSequence
from astropy.nddata import StdDevUncertainty
from sunpy.visualization.colormaps import color_tables as ct
from EMToolKit.algorithms.simple_reg_dem_wrapper import autoloading_simple_reg_dem_wrapper

def multi_inst_simple_dem_wrapper(datasequence, wrapargs=None, do_plot=False, data_dir="../data/multitest", recalc_simple=False):
    """
    Perform simple DEM analysis using multiple instruments.

    Parameters:
        datasequence (NDCubeSequence): Sequence of data cubes for multiple instruments.
        wrapargs (dict, optional): Arguments for the DEM wrapper. Defaults to None.
        do_plot (bool, optional): Whether to plot the reprojected data immediately. Defaults to True.
        data_dir (str, optional): Directory for data storage. Defaults to "../data/multitest".
        recalc_simple (bool, optional): Whether to recalculate simple DEM. Defaults to False.

    Returns:
        NDCubeSequence: Resulting DEM analysis.
    """
    if wrapargs is None:
        wrapargs = {}

    downprojected_sequence, nan_mask, coarsest_cube = reproject_with_uncertainties(datasequence)

    if do_plot:
        plot_reprojected(downprojected_sequence, nan_mask, coarsest_cube)

    wrapargs["prepend"] = "multi_inst_"

    return autoloading_simple_reg_dem_wrapper(downprojected_sequence, data_dir, wrapargs=wrapargs, recalc_simple=recalc_simple)

def find_coarsest_cube(datasequence):
    """
    Find the cube with the largest pixel scale in the sequence.

    Parameters:
        datasequence (NDCubeSequence): Sequence of data cubes for multiple instruments.

    Returns:
        NDCube: The cube with the largest pixel scale.
    """
    max_scale = 0
    coarsest_cube = None
    for seq in datasequence:
        pixel_scale = seq.meta['CDELT1'] * seq.meta['CDELT2']
        if pixel_scale >= max_scale:
            max_scale = pixel_scale
            coarsest_cube = seq
    if coarsest_cube is None:
        raise ValueError("No coarsest cube found in the sequence")
    return coarsest_cube

def reproject_cubes(fine_sequence, coarsest_cube):
    """
    Reproject fine maps to the coarse map shape.

    Parameters:
        fine_sequence (NDCubeSequence): Sequence of fine data cubes.
        coarsest_cube (NDCube): The coarsest data cube.

    Returns:
        NDCubeSequence: Sequence of reprojected data cubes.
    """
    reprojected_cubes = [mp.reproject_to(coarsest_cube.wcs) for mp in fine_sequence]
    return NDCubeSequence(reprojected_cubes, meta=fine_sequence.meta)

def scale_uncertainties(fine_sequence, downprojected_sequence, coarsest_cube):
    """
    Scale and reproject uncertainties.

    Parameters:
        fine_sequence (NDCubeSequence): Sequence of fine data cubes.
        downprojected_sequence (NDCubeSequence): Sequence of reprojected data cubes.
        coarsest_cube (NDCube): The coarsest data cube.

    Returns:
        list: List of reprojected uncertainties.
    """
    orig_area = fine_sequence[0].data.shape[0] * fine_sequence[0].data.shape[1]
    new_area = downprojected_sequence[0].data.shape[0] * downprojected_sequence[0].data.shape[1]
    area_ratio_rt = np.sqrt(new_area / orig_area)

    uncertainties = [
        StdDevUncertainty(
            NDCube(
                area_ratio_rt * mp.uncertainty.array.data, mp.wcs, meta=mp.meta
            ).reproject_to(coarsest_cube.wcs).data
        ) for mp in fine_sequence
    ]
    return uncertainties

def combine_sequences(downprojected_sequence, coarsest_cube, uncertainties, nan_level):
    """
    Combine reprojected fine maps with the coarse map and apply NaN mask.

    Parameters:
        downprojected_sequence (NDCubeSequence): Sequence of reprojected data cubes.
        coarsest_cube (NDCube): The coarsest data cube.
        uncertainties (list): List of reprojected uncertainties.
        nan_level (float): Threshold level for NaN masking.

    Returns:
        NDCubeSequence: Combined sequence of data cubes.
        numpy.ndarray: NaN mask.
    """
    nan_mask = np.where(np.isnan(uncertainties[0].array), np.nan, 1)
    for cube in downprojected_sequence:
        nan_mask *= np.where(np.isnan(cube.data) | (cube.data < nan_level), np.nan, 1)

    combined_cubes = [NDCube(cube.data * nan_mask, meta=cube.meta, wcs=cube.wcs, uncertainty=cube.uncertainty) for cube in downprojected_sequence]
    combined_cubes.append(NDCube(coarsest_cube.data * nan_mask, meta=coarsest_cube.meta, wcs=coarsest_cube.wcs))

    return NDCubeSequence(combined_cubes, common_axis=0, meta=downprojected_sequence.meta), nan_mask

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

    coarsest_cube = find_coarsest_cube(datasequence)
    fine_sequence = NDCubeSequence([mp for mp in datasequence if mp is not coarsest_cube], meta=datasequence.meta)

    downprojected_sequence = reproject_cubes(fine_sequence, coarsest_cube)
    uncertainties = scale_uncertainties(fine_sequence, downprojected_sequence, coarsest_cube)

    for i in range(len(downprojected_sequence)):
        downprojected_sequence[i].uncertainty = uncertainties[i]

    combined_sequence, nan_mask = combine_sequences(downprojected_sequence, coarsest_cube, uncertainties, nan_level)

    print("Reprojected successfully!")

    return combined_sequence, nan_mask, coarsest_cube

def plot_reprojected(downprojected_sequence, nan_mask, coarsest_cube):
    """
    Plot the reprojected maps and their uncertainties.

    Parameters:
        downprojected_sequence (NDCubeSequence): Sequence of reprojected data cubes.
        nan_mask (numpy.ndarray): Mask for NaN values.
        coarsest_cube (NDCube): Coarsest data cube for comparison.
    """
    for aia_reproj_map in downprojected_sequence:
        if aia_reproj_map is coarsest_cube:
            continue

        fig = plt.figure(figsize=(12, 6))

        # Apply NaN mask to the data
        coarsest_data_masked = coarsest_cube.data * nan_mask
        aia_data_masked = aia_reproj_map.data * nan_mask

        # Create axes with WCS projection
        ax1 = fig.add_subplot(1, 2, 1, projection=aia_reproj_map.wcs)
        ax2 = fig.add_subplot(1, 2, 2, projection=coarsest_cube.wcs, sharex=ax1, sharey=ax1)

        ax1.set_facecolor("grey")
        ax2.set_facecolor("grey")

        ax1.grid(True)
        ax2.grid(True)

        # Plot the coarsest cube data
        ax2.imshow(coarsest_data_masked, cmap='gray', origin="lower")

        # Get wavelength for the AIA map
        wave = aia_reproj_map.meta.get("wavelnth", None)
        if wave is None:
            continue
        else:
            wave *= u.angstrom

        # Overlay the AIA reprojected data
        ax1.imshow(np.sqrt(aia_data_masked), alpha=0.75, cmap=ct.aia_color_table(wave))

        ax1.set_title(f'AIA {wave}; Shape: {aia_reproj_map.data.shape}')
        ax2.set_title(f'XRT Shape: {coarsest_cube.data.shape}')

        plt.tight_layout()
        plt.show()
