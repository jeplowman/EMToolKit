import copy, numpy as np
from sunpy.map import Map
from ndcube import NDCube, NDCubeSequence, NDCollection
from astropy.coordinates import SkyCoord
from astropy.nddata import StdDevUncertainty

# The code here is only a prototype/placeholder!

# Given a set of XRT SunPy Maps, return the appropriate arguments for use
# as an EMToolKit data sequence -- the selection of maps appropriate for
# DEMs, corresponding errors, temperature response
# functions and corresponding (log) temperature arrays
def xrt_wrapper(maps_in, temperature_array):
	[maps,logts,tresps,errs] = [[],[],[],[]]
	for i in range(0,len(maps_in)):
		current_map = copy.deepcopy(maps_in[i])#.rotate(order=3)
		current_map.meta['channel'] = current_map.meta['ec_fw1_']+'_'+current_map.meta['ec_fw2_']
		if(not('detector' in current_map.meta)): current_map.meta['detector'] = 'XRT'
		[logt,tresp] = xrt_temperature_response(current_map, temperature_array=temperature_array)
		if(len(tresp) == len(logt)):
			maps.append(current_map)
			errs.append(StdDevUncertainty(estimate_xrt_error(current_map)))
			logts.append(logt)
			tresps.append(tresp)
	return maps,errs,logts,tresps

def estimate_xrt_error(map_in):
	channel = map_in.meta['ec_fw1_']+'_'+map_in.meta['ec_fw2_']
	if(map_in.meta['ec_fw1_'] != 'Open' and map_in.meta['ec_fw2_'] != 'Ti_poly'):
		print("Estimate XRT error is a prototype and only empirically estimates Open Ti-poly uncertainties")
	return ((np.abs(map_in.data.astype(np.float32))*25 + 25**2)**0.5)





def xrt_temperature_response(map_in, temperature_array):
	import xrtpy
	import matplotlib.pyplot as plt
	import numpy as np
	# The real temperature response functions are found here:
	# https://xrtpy.readthedocs.io/en/stable/notebooks/computing_functions/temperature_response.html

	# Available channels: "Al-mesh", "Al-poly", "C-poly", "Ti-poly", "Be-thin", "Be-med", "Al-med", "Al-thick", "Be-thick" , "Al-poly/Al-mesh", "Al-poly/Ti-poly", "Al-poly/Al-thick", "Al-poly/Be-thick" , "C-poly/Ti-poly"

	channel = map_in.meta['ec_fw1_']+'/'+map_in.meta['ec_fw2_']
	channel = channel.replace("Open/", "")
	print("MAP FILTER IS : ", channel)
	filter = channel
	date_time = "2007-09-22T21:59:59"
	Temperature_Response_Fundamental = xrtpy.response.TemperatureResponseFundamental(
		filter, date_time, abundance_model="Photospheric"
	)
	temperature_response = Temperature_Response_Fundamental.temperature_response()
	CHIANTI_temperature = Temperature_Response_Fundamental.CHIANTI_temperature
	log_CHIANTI_temperature = np.log10(CHIANTI_temperature.value)
	logt, tresp = interp1d_logt(log_CHIANTI_temperature, temperature_response.value, temperature_array)
	plt.plot(logt, tresp, 'b')
	plt.show()
	print(len(logt), len(tresp))
	print(type(logt), type(tresp))
	return logt, tresp


def xrt_temperature_response_old(map_in, step_size=0.2):
	import xrtpy
	import matplotlib.pyplot as plt
	import numpy as np
	# The real temperature response functions are found here:
	# https://xrtpy.readthedocs.io/en/stable/notebooks/computing_functions/temperature_response.html

	channel = map_in.meta['ec_fw1_']+'_'+map_in.meta['ec_fw2_']
	refchannels = np.array(['Open_Ti_poly'])

	logt_table = np.array([5.5, 6.0, 6.2, 6.4, 6.5, 6.6, 6.8, 6.95, 7.1, 7.2, 7.35, 7.5])
	# This preliminary temperature response for the XRT open/Ti-poly channel is from visually reading off
	# The figure at http://solar.physics.montana.edu/takeda/xrt_response/comp_xre2_Tipol.png
	tresp_table = np.array([[0.0], [0.5], [1.0], [2.0], [4.0], [6.0], [15.], [16.5], [15.], [7.0], [4.5], [3.5]])*1.0e-26

	logt, tresp = interpolate_table(logt_table, tresp_table, step_size)
	resp_out = tresp[:,np.where(refchannels == channel)].flatten()

	plt.plot(logt, resp_out, 'r')
	plt.show()

	print(len(logt), len(resp_out))
	print(type(logt), type(resp_out))
	return logt, resp_out


def interp1d_logt(logt, tresp, temperature_array):
	new_tresp = np.interp(temperature_array, logt, tresp)
	return temperature_array, new_tresp


def interpolate_table(logt, table, step_size):
    new_logt = np.arange(np.min(logt), np.max(logt) + step_size, step_size)

    interpolated_table = np.zeros((len(new_logt), table.shape[1]))

    for i in range(table.shape[1]):
        interpolated_table[:, i] = np.interp(new_logt, logt, table[:, i])

    return new_logt, interpolated_table


import os
import numpy as np
import astropy.units as u
from sunpy.net import Fido, attrs as a
from sunpy.time import TimeRange
from astropy.time import Time
from EMToolKit.util import list_fits_files


def download_xrt_data(base_path, date, redownload=False):
	folder_name = date.replace("/", "_").replace(" ", "_").replace(":", "_")
	xrt_data_dir = os.path.join(base_path, ".data", folder_name)

	if not os.path.exists(xrt_data_dir):
		os.makedirs(xrt_data_dir)

	paths = list_fits_files(xrt_data_dir, 'xrt')

	print(f"Found {len(paths)} xrt images on disk.")

	if len(paths) == 0 or redownload:
		try:
			print(f"Searching for XRT images from {date}...")
			qry = Fido.search(a.Time(TimeRange(date, 4 * u.s)), a.Instrument('XRT'))
			print(f"Downloading XRT images...")
			Fido.fetch(qry, path=xrt_data_dir, max_conn=3)
		except ConnectionError:
			print(f"No internet connection, continuing without xrt")
		paths = list_fits_files(xrt_data_dir, 'xrt')
	return paths, xrt_data_dir