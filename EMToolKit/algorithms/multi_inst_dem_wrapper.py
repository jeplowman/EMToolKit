
import matplotlib.pyplot as plt
import astropy.units as u
import os, numpy as np

from ndcube import NDCube, NDCubeSequence, NDCollection
from sunpy.visualization.colormaps.color_tables import aia_color_table, xrt_color_table
from astropy.nddata import StdDevUncertainty

# datasequence can behave in a list-like fashion with
# n elements that behave as sunpy maps that have the following
# 	data
# 	uncertainty.array (same shape as data)
# 	wcs
#	observer_coordinate
# meta including temperature response, log temperature array
# and exposure time
# The spatial response of each pixel must be either supplied
# or estimated based on the wcs. For instruments that do not
# have spatially localized detector elements (e.g., RHESSI),
# this will presumably have to be supplied.
# For instruments where the spatial and temperature response
# are not separable (e.g., overlappographs), this will need
# to be indicated in the meta, and the responses will
# need to be supplied, methods for computed them provided, or
# a standard way of defining them must be developed.
# The current baseline in EMToolKit assumes localized detector
# elements and only provides wcs for spatial information.
# We should begin by implementing a backward compatible layer
# for that.

def multi_inst_simple_dem_wrapper(datasequence, wrapargs={}, doPlot=False, dat_dir="../data/multitest", recalc_simple=False):

	downprojected_sequence, nan_mask, coarsest_cube= _reproject_with_uncertainties(datasequence)

	if doPlot:
		plot_reprojected(downprojected_sequence, nan_mask, coarsest_cube)

	print("Performing simple DEM...")
	wrapargs["prepend"]="multi_inst_"
	from EMToolKit.algorithms.simple_reg_dem_wrapper import autoloading_simple_reg_dem_wrapper
	return autoloading_simple_reg_dem_wrapper(datasequence, dat_dir, wrapargs=wrapargs, recalc_simple=recalc_simple)


def _reproject_with_uncertainties(datasequence, nan_level=-50):

		print("Reprojecting all maps to smallest map...", end="")
		# Find the map with the largest pixel scale
		sz_min = 0
		for seq in datasequence:
			sz = seq.meta['CDELT1']* seq.meta['CDELT2']
			if sz >= sz_min:
				sz_min = sz
				coarsest_cube = seq

		# Pull out the fine Sequence
		fine_sequence = NDCubeSequence([mp for mp in datasequence if mp is not coarsest_cube], meta=datasequence.meta)

		# Reproject the fine maps to the coarse map shape
		downprojected_sequence = NDCubeSequence([mp.reproject_to(coarsest_cube.wcs) for mp in fine_sequence])

		# Compute factor to scale the uncertainties by the area ratio
		orig_area = fine_sequence[0].data.shape[0]*fine_sequence[0].data.shape[1]
		new_area = downprojected_sequence[0].data.shape[0]*downprojected_sequence[0].data.shape[1]
		area_ratio_rt = np.sqrt(new_area/orig_area)

		# Reproject and scale the uncertainties
		uncertainties = [StdDevUncertainty(
			NDCube(area_ratio_rt*mp.uncertainty.array.data,  mp.wcs,  meta=mp.meta).
			reproject_to(coarsest_cube.wcs).data
											) for mp in fine_sequence]
		for i in range(len(downprojected_sequence)):
			downprojected_sequence[i].uncertainty = uncertainties[i]

		# Combine the reprojected fine maps with the coarse map
		nan_mask = np.where(np.isnan(uncertainties[0].array), np.nan, 1)
		for cube in downprojected_sequence:
			nan_mask *= np.where(cube.data<nan_level, np.nan, 1)
		full_list = [x * nan_mask for x in downprojected_sequence]
		full_list.append(coarsest_cube * nan_mask)
		datasequence = NDCubeSequence(full_list, common_axis=0)

		print("Reprojected successfully!")

		return datasequence, nan_mask, coarsest_cube

def plot_reprojected(downprojected_sequence, nan_mask, coarsest_cube):
	# Plot the reprojected maps
	for aia_reproj_map in downprojected_sequence:
		fig = plt.figure(figsize=(12, 6))
		coarsest_cube *= nan_mask
		aia_reproj_map *= nan_mask
		ax1 = fig.add_subplot(1, 2, 1, projection=coarsest_cube)
		ax2 = fig.add_subplot(1, 2, 2, projection=coarsest_cube, sharex=ax1, sharey=ax1)
		ax1.set_facecolor("grey")
		ax2.set_facecolor("grey")
		ax1.imshow(coarsest_cube.data, cmap=xrt_color_table(), origin="lower")
		wave = aia_reproj_map.meta.get("wavelnth", 0)*u.angstrom
		# aia_reproj_map.plot( axes=ax1, alpha=0.75, cmap=aia_color_table(wave))
		ax1.imshow(np.sqrt(aia_reproj_map.data), alpha=0.75, cmap=aia_color_table(wave))
		ax2.imshow(aia_reproj_map.uncertainty.array, cmap="plasma")
		ax1.set_title(f'AIA {wave} overlaid on XRT\nshape: {aia_reproj_map.data.shape}')
		ax2.set_title(f'Uncertainty\nshape: {aia_reproj_map.uncertainty.array.shape}')
		plt.tight_layout()
		plt.show()