
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
def multi_inst_simple_dem_wrapper(datasequence, wrapargs={}, doPlot=False, dat_dir="../data/multitest"):

	print("Reprojecting all maps to smallest map...", end="")
	# Find the smallest map, which should be the XRT map
	sz_min = 1e9
	for seq in datasequence:
		sz = seq.data.shape[0] * seq.data.shape[1]
		if sz <= sz_min:
			sz_min = sz
			xrt_sequence = seq

	# Pull out the AIA Sequence
	aia_sequence = NDCubeSequence([mp for mp in datasequence if mp is not xrt_sequence], meta=datasequence.meta)

	# Reproject the AIA maps to the XRT map shape
	aia_reproj_seqs = NDCubeSequence([mp.reproject_to(xrt_sequence.wcs) for mp in aia_sequence])

	# Compute factor to scale the uncertainties by the area ratio
	orig_shape = aia_sequence[0].data.shape
	orig_area = orig_shape[0]*orig_shape[1]
	new_shape = aia_reproj_seqs[0].data.shape
	new_area = new_shape[0]*new_shape[1]
	area_ratio_rt = np.sqrt(new_area/orig_area)

	# Reproject and scale the uncertainties
	uncertainties = [StdDevUncertainty(NDCube(
		area_ratio_rt*mp.uncertainty.array.data, mp.wcs, meta=mp.meta).
		reproject_to(xrt_sequence.wcs).data) for mp in aia_sequence]

	for i in range(len(aia_reproj_seqs)):
		aia_reproj_seqs[i].uncertainty = uncertainties[i]

	# Combine the reprojected AIA maps with the XRT map
	aia_list = [x for x in aia_reproj_seqs]
	aia_list.append(xrt_sequence)
	new_sequence = NDCubeSequence(aia_list, common_axis=0)
	datasequence = new_sequence

	if doPlot:
		# Plot the reprojected maps
		for aia_reproj_map in aia_reproj_seqs:
			print(type(aia_reproj_map))
			fig = plt.figure(figsize=(12, 6))
			ax1 = fig.add_subplot(1, 2, 1, projection=xrt_sequence)
			ax2 = fig.add_subplot(1, 2, 2, projection=xrt_sequence)
			xrt_sequence.plot(axes=ax1, cmap=xrt_color_table())
			wave = aia_reproj_map.meta["wavelnth"]*u.angstrom
			aia_reproj_map.plot( axes=ax1, alpha=0.5, cmap=aia_color_table(wave))
			ax2.imshow(aia_reproj_map.uncertainty.array, cmap="plasma")
			ax1.set_title(f'AIA {wave} overlaid on XRT\nshape: {aia_reproj_map.data.shape}')
			ax2.set_title(f'Uncertainty\nshape: {aia_reproj_map.uncertainty.array.shape}')
			plt.tight_layout()
			plt.show()

	print("Success!")
	print("Performing simple DEM...")
	from EMToolKit.algorithms.simple_reg_dem_wrapper import autoloading_simple_reg_dem_wrapper
	return autoloading_simple_reg_dem_wrapper(datasequence, dat_dir, wrapargs=wrapargs)

