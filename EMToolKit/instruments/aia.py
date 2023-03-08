import copy, numpy as np
from sunpy.map import Map
from ndcube import NDCube, NDCubeSequence, NDCollection
from astropy.coordinates import SkyCoord
from astropy.nddata import StdDevUncertainty

# H. Farooki: Temporarily set missing=0 because it was giving me an error otherwise

# Load AIA data from a set of paths and return the appropriate 
# arguments for use in an EMToolKit DataSequence -- a list
# of SunPy maps, corresponding errors, the log temperature
# axes for the temperature response functions and the temperature
# response functions themselves:
def load_from_paths(paths,xl=None,yl=None,dx=None,dy=None,refindex=0):
	refmap = Map(paths[refindex]).rotate(order=3, missing=0)
	nocrop = (xl is None or yl is None or dx is None or dy is None)
	if(nocrop == False):
		blc=SkyCoord(xl,yl,frame=refmap.coordinate_frame)
		trc=SkyCoord(xl+dx,yl+dy,frame=refmap.coordinate_frame)

	maps=[]
	for i in range(0,len(paths)):
		maps.append(Map(paths[i]).rotate(order=3, missing=0))
		if(nocrop==False): maps[i] = maps[i].submap(blc,top_right=trc)
	return maps

# Given a set up AIA SunPy Maps, return the appropriate arguments for use 
# as an EMToolKit data sequence -- the selection of maps appropriate for
# DEMs (EUV not including 304), corresponding errors, temperature response
# functions and corresponding (log) temperature arrays
def aia_wrapper(maps_in):
	[maps,logts,tresps,errs] = [[],[],[],[]]
	for i in range(0,len(maps_in)):
		current_map = copy.deepcopy(maps_in[i]).rotate(order=3, missing=0)	
		if(not('detector' in current_map.meta)): current_map.meta['detector'] = 'AIA'
		[logt,tresp] = aia_temperature_response(current_map)
		if(len(tresp) == len(logt)):
			maps.append(current_map)
			errs.append(StdDevUncertainty(estimate_aia_error(current_map)))
			logts.append(logt)
			tresps.append(tresp)
	return maps,errs,logts,tresps

# These are placeholder methods pending temperature response functions
# being added to AIApy. error estimates are included for the initial 
# version as well, although I believe they're already in the initial version
def estimate_aia_error(map_in):
	channel = map_in.meta['detector']+map_in.meta['wave_str']
	refchannels=np.array(['AIA94_THIN', 'AIA131_THIN', 'AIA171_THIN', 'AIA193_THIN', 'AIA211_THIN', 'AIA304_THIN', 'AIA335_THIN'])
	refg = np.array([2.128,1.523,1.168,1.024,0.946,0.658,0.596])
	refn = np.array([1.14,1.18,1.15,1.2,1.2,1.14,1.18])
	sigmas = np.zeros(map_in.data.shape)
	dnpp = refg[np.where(refchannels == channel)]
	rdn = refn[np.where(refchannels == channel)]
	return np.sqrt(np.clip(map_in.data*dnpp,0.0,None) + rdn**2)
	
def aia_temperature_response(map_in):

	import numpy as np
	channel = map_in.meta['detector']+map_in.meta['wave_str']
	refchannels=np.array(['AIA94_THIN', 'AIA131_THIN', 'AIA171_THIN', 'AIA193_THIN', 'AIA211_THIN', 'AIA335_THIN'])

	tresp_table = np.array([
		[7.40163203e-29, 1.46072515e-26, 5.34888608e-26, 8.50682605e-26, 2.03073186e-26, 3.35720499e-27],
		[1.40367160e-28, 2.45282476e-26, 5.69953432e-26, 7.63243171e-26, 1.69626905e-26, 2.63732660e-27],
		[2.34619332e-28, 3.55303726e-26, 8.99381587e-26, 8.00857294e-26, 1.61667405e-26, 2.65194161e-27],
		[3.50314745e-28, 4.59075911e-26, 1.69634100e-25, 8.95460085e-26, 1.65561573e-26, 2.97594545e-27],
		[4.78271991e-28, 5.36503515e-26, 3.12477095e-25, 9.96536550e-26, 1.75607919e-26, 3.35107900e-27],
		[6.12070853e-28, 5.68444094e-26, 5.24806168e-25, 1.08520987e-25, 1.91638109e-26, 3.66488663e-27],
		[7.55587098e-28, 5.47222296e-26, 7.88759175e-25, 1.16018001e-25, 2.13919425e-26, 3.86907993e-27],
		[9.17523861e-28, 4.81119761e-26, 1.05392610e-24, 1.23253858e-25, 2.41030560e-26, 3.98249755e-27],
		[1.08989852e-27, 3.85959059e-26, 1.24465730e-24, 1.33290059e-25, 2.69021271e-26, 4.00495098e-27],
		[1.24045630e-27, 2.80383406e-26, 1.28593896e-24, 1.53051294e-25, 2.93365234e-26, 3.88113016e-27],
		[1.32502623e-27, 1.83977650e-26, 1.14316216e-24, 1.96447388e-25, 3.18854130e-26, 3.62558008e-27],
		[1.31341529e-27, 1.11182849e-26, 8.56234347e-25, 2.81167898e-25, 3.70341037e-26, 3.33826036e-27],
		[1.20609026e-27, 6.50748885e-27, 5.28341859e-25, 4.06492781e-25, 4.99880175e-26, 3.04882976e-27],
		[1.04215049e-27, 3.96843481e-27, 2.62971839e-25, 5.18462720e-25, 7.85106324e-26, 2.71070234e-27],
		[9.10286530e-28, 2.61876319e-27, 1.06118152e-25, 5.21700296e-25, 1.24919334e-25, 2.50715036e-27],
		[8.66408740e-28, 1.85525324e-27, 4.07695483e-26, 3.80438138e-25, 1.66057469e-25, 2.81207972e-27],
		[7.91289650e-28, 1.39717024e-27, 2.52503209e-26, 1.89390887e-25, 1.57937960e-25, 3.70484370e-27],
		[5.79545741e-28, 1.11504283e-27, 2.46082719e-26, 6.76946014e-26, 1.03164467e-25, 4.56853237e-27],
		[3.76980809e-28, 9.38169611e-28, 1.99231167e-26, 2.29005173e-26, 5.22396021e-26, 4.79107749e-27],
		[3.07180587e-28, 8.24801234e-28, 1.03293005e-26, 1.04739857e-26, 2.40915758e-26, 4.45474563e-27],
		[3.80396260e-28, 7.43331919e-28, 3.95012736e-27, 7.10153621e-27, 1.12964568e-26, 3.86954529e-27],
		[5.90602903e-28, 6.74537063e-28, 1.80531856e-27, 6.11725282e-27, 5.94188041e-27, 3.22466361e-27],
		[9.44297921e-28, 6.14495760e-28, 1.34240926e-27, 5.72716914e-27, 3.88586017e-27, 2.60315761e-27],
		[1.44620154e-27, 5.70805277e-28, 1.22599436e-27, 5.81269453e-27, 3.15819472e-27, 2.03714253e-27],
		[2.07453706e-27, 5.61219786e-28, 1.05241421e-27, 6.32309411e-27, 2.71019768e-27, 1.53754765e-27],
		[2.75538782e-27, 6.31981777e-28, 8.14794974e-28, 6.39019697e-27, 2.11642133e-27, 1.11256462e-27],
		[3.34721420e-27, 9.19747307e-28, 6.22671197e-28, 5.43856741e-27, 1.48862847e-27, 7.71472387e-28],
		[3.65940611e-27, 1.76795732e-27, 5.15934126e-28, 4.02432562e-27, 1.02304623e-27, 5.25924670e-28],
		[3.52350628e-27, 3.77985446e-27, 4.75354344e-28, 2.82301679e-27, 7.31784516e-28, 3.94053508e-28],
		[2.91029046e-27, 7.43166191e-27, 4.74903304e-28, 2.11704426e-27, 5.53478904e-28, 3.84310717e-28],
		[2.00104132e-27, 1.19785603e-26, 4.78687900e-28, 2.23337165e-27, 4.46856485e-28, 4.58924041e-28],
		[1.11792532e-27, 1.48234676e-26, 4.55829373e-28, 3.97858677e-27, 3.92863480e-28, 5.15941951e-28],
		[5.06298440e-28, 1.36673114e-26, 4.07150741e-28, 7.95264076e-27, 3.72931301e-28, 4.64395852e-28],
		[1.98518377e-28, 9.61047146e-27, 3.57245583e-28, 1.28401071e-26, 3.65339284e-28, 3.28256715e-28],
		[8.04488438e-29, 5.61209353e-27, 3.20550268e-28, 1.62347002e-26, 3.55630868e-28, 1.96954455e-28],
		[4.22436839e-29, 3.04779780e-27, 2.95165454e-28, 1.72337033e-26, 3.40535153e-28, 1.12675591e-28],
		[3.04428119e-29, 1.69378976e-27, 2.76040754e-28, 1.64931273e-26, 3.22463539e-28, 6.77096153e-29],
		[2.64974316e-29, 1.02113491e-27, 2.60329069e-28, 1.49059780e-26, 3.03925227e-28, 4.49240089e-29],
		[2.48278975e-29, 6.82223774e-28, 2.46810279e-28, 1.30615389e-26, 2.86214839e-28, 3.30827404e-29],
		[2.38426184e-29, 5.02099099e-28, 2.34850846e-28, 1.12339098e-26, 2.69680531e-28, 2.65114499e-29],
		[2.30756823e-29, 3.99377760e-28, 2.24072357e-28, 9.54159826e-27, 2.54348518e-28, 2.25550937e-29]
	])

	logt = np.array([5.5 , 5.55, 5.6 , 5.65, 5.7 , 5.75, 5.8 , 5.85, 5.9 , 5.95, 6., 
					6.05, 6.1 , 6.15, 6.2 , 6.25, 6.3 , 6.35, 6.4 , 6.45, 6.5 , 
					6.55, 6.6 , 6.65, 6.7 , 6.75, 6.8 , 6.85, 6.9 , 6.95, 7.  , 
					7.05, 7.1 , 7.15, 7.2 , 7.25, 7.3 , 7.35, 7.4 , 7.45, 7.5 ])
					
	return logt,tresp_table[:,np.where(refchannels == channel)].flatten()
