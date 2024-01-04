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
def xrt_wrapper(maps_in):
	[maps,logts,tresps,errs] = [[],[],[],[]]
	for i in range(0,len(maps_in)):
		current_map = copy.deepcopy(maps_in[i])#.rotate(order=3)	
		current_map.meta['channel'] = current_map.meta['ec_fw1_']+'_'+current_map.meta['ec_fw2_']
		if(not('detector' in current_map.meta)): current_map.meta['detector'] = 'XRT'
		[logt,tresp] = xrt_temperature_response(current_map)
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

def xrt_temperature_response(map_in):
	channel = map_in.meta['ec_fw1_']+'_'+map_in.meta['ec_fw2_']
	refchannels = np.array(['Open_Ti_poly'])
	logt = np.array([5.5, 6.0, 6.2, 6.4, 6.5, 6.6, 6.8, 6.95, 7.1, 7.2, 7.35, 7.5])
	# This preliminary temperature response for the XRT open/Ti-poly channel is from visually reading off
	# The figure at http://solar.physics.montana.edu/takeda/xrt_response/comp_xre2_Tipol.png
	tresp_table = np.array([[0.0], [0.5], [1.0], [2.0], [4.0], [6.0], [15.], [16.5], [15.], [7.0], [4.5], [3.5]])*1.0e-26
	
	return logt,tresp_table[:,np.where(refchannels == channel)].flatten()
	
