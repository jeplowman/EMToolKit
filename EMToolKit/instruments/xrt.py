import copy, numpy as np
from sunpy.map import Map
from ndcube import NDCube, NDCubeSequence, NDCollection
from astropy.coordinates import SkyCoord
from astropy.nddata import StdDevUncertainty
import xrtpy
import matplotlib.pyplot as plt
import os
import numpy as np
import astropy.units as u
from sunpy.net import Fido, attrs as a
from sunpy.time import TimeRange
from astropy.time import Time
from EMToolKit.util import list_fits_files

# Given a set of XRT SunPy Maps, return the appropriate arguments for use
# as an EMToolKit data sequence -- the selection of maps appropriate for
# DEMs, corresponding errors, temperature response
# functions and corresponding (log) temperature arrays
def xrt_wrapper(maps_in, temperature_array):
	"""
	Processes a list of XRT maps and computes their temperature responses.
	This function returns the modified maps along with their associated uncertainties and temperature response data.

	Args:
		maps_in (list): A list of input maps to be processed.
		temperature_array (array-like): An array of temperatures used for calculating the temperature response.

	Returns:
		tuple: A tuple containing the processed maps, their uncertainties, logarithmic temperature values, and temperature responses.

	Raises:
		ValueError: If the input maps are not compatible with the temperature array.
	"""

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


def xrt_temperature_response(map_in, temperature_array, *, do_plot=False):
	"""
	This Python function calculates the temperature response for a given filter channel using xrtpy
	library and plots the response curve.

	Args:
		map_in: The sunpy map to get the channel information from
		temperature_array: The `temperature_array` parameter in the `xrt_temperature_response` function is
		an array that contains the temperature values for which you want to calculate the temperature
		response. This array will be used as input to interpolate the temperature response values obtained
		from the `Temperature_Response_Fundamental` object.

	Returns:
		logt: the logarithm of CHIANTI temperature values
  		tresp: the corresponding temperature response values
	"""

	# The real temperature response functions are found here:
	# https://xrtpy.readthedocs.io/en/stable/notebooks/computing_functions/temperature_response.html

	# Available channels: "Al-mesh", "Al-poly", "C-poly", "Ti-poly", "Be-thin", "Be-med", "Al-med", "Al-thick",
 	# "Be-thick" , "Al-poly/Al-mesh", "Al-poly/Ti-poly", "Al-poly/Al-thick", "Al-poly/Be-thick" , "C-poly/Ti-poly"

	channel = map_in.meta['ec_fw1_']+'/'+map_in.meta['ec_fw2_']
	channel = channel.replace("Open/", "")
	filt = channel
	date_time = "2007-09-22T21:59:59"
	Temperature_Response_Fundamental = xrtpy.response.TemperatureResponseFundamental(
		filt, date_time, abundance_model="Photospheric"
	)
	temperature_response = Temperature_Response_Fundamental.temperature_response()
	CHIANTI_temperature = Temperature_Response_Fundamental.CHIANTI_temperature
	log_CHIANTI_temperature = np.log10(CHIANTI_temperature.value)
	logt, tresp = interp1d_logt(log_CHIANTI_temperature, temperature_response.value, temperature_array)
	if do_plot:
		plt.plot(logt,tresp)
		plt.title(f"XRT MAP FILTER IS : {channel}")
		plt.xlabel("Temperature")
		plt.show()
	return logt, tresp


def interp1d_logt(logt, tresp, temperature_array):
	"""
	Interpolates the temperature response values for a given set of logarithmic temperature values.
	This function maps the temperature response to the specified temperature array using linear interpolation.

	Args:
		logt (array-like): An array of logarithmic temperature values.
		tresp (array-like): An array of temperature response values corresponding to the logarithmic temperatures.
		temperature_array (array-like): The array of temperatures for which the response is to be interpolated.

	Returns:
		tuple: A tuple containing the input temperature array and the interpolated temperature response values.
	"""

	new_tresp = np.interp(temperature_array, logt, tresp)
	return temperature_array, new_tresp



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
			print("Downloading XRT images...")
			Fido.fetch(qry, path=xrt_data_dir, max_conn=3)
		except ConnectionError:
			print("No internet connection, continuing without xrt")
		paths = list_fits_files(xrt_data_dir, 'xrt')
	return paths, xrt_data_dir