import copy
import numpy as np
from sunpy.map import Map
from ndcube import NDCube, NDCubeSequence, NDCollection
from astropy.coordinates import SkyCoord
from astropy.nddata import StdDevUncertainty
from astropy.nddata import UnknownUncertainty
from emtoolkit.schemas.basic_schemas import basic_detector, basic_source
import astropy

def em_data(maps,errs,logts,tresps,channels=None):
	cubes = []
	for i in range(0,len(maps)):
		mapi = copy.copy(maps[i])
		if('CHANNEL' not in mapi.meta):
			if(channels is None): mapi.meta['CHANNEL'] = mapi.meta['DETECTOR']+mapi.meta['WAVE_STR']
			else: mapi.meta['CHANNEL'] = channels[i]
		mapi.meta['LOGT'] = logts[i]
		mapi.meta['TRESP'] = tresps[i]
		mapi.meta['SHAPE'] = np.array(mapi.data.shape)
		if('SCHEMA' not in mapi.meta):
			mapi.meta['SCHEMA'] = basic_detector(mapi.meta)
		cubes.append(NDCube(mapi,uncertainty=errs[i]))
	return NDCubeSequence(cubes)

def dem_model(coeffs,logts,bases,coord_info,algorithm,wrapper,meta=None,wrapargs=None):

	nd = len(coeffs)
	dem_sequence = []
	if(isinstance(coord_info,astropy.wcs.wcs.WCS)): wcs = coord_info
	else: # Check to see if the coordinate info input is a dict:
		meta = coord_info
		if(isinstance(coord_info,dict)==False):
			print('Warning in emtoolkit dem_model: coord_info is not a wcs or dict')
		wcs = meta.get('wcs',None)

	if(meta is None):
		if(wcs is not None): meta = dict(wcs.to_header())
		else: print('Warning in emtoolkit dem_model: need wcs or image meta')
	if('LOGT' not in meta): meta['LOGT'] = logts[0]
	if('SCHEMA' not in meta):
		meta['SCHEMA'] = basic_source([Map(coeffs[0],meta)])
	for i in range(0,nd):
		logt0 = np.median(logts[i][np.where(bases[i]==np.max(bases[i]))])
		metai = {'BASIS':bases[i],'LOGT':logts[i],'LOGT0':logt0,'ALGORITHM':algorithm,'WRAPPER':wrapper,'WRAPARGS':wrapargs}
		dem_sequence.append(NDCube(coeffs[i],wcs=wcs,meta={**meta,**metai}))
	return NDCubeSequence(dem_sequence,meta={'ALGORITHM':algorithm,'WRAPPER':wrapper,'WRAPARGS':wrapargs})



class em_collection:

	def __init__(self,datasequence):
		self.collection = NDCollection([("data",datasequence),("models",[])])

	def data(self): return self.collection['data']

	def add_model(self,modelsequence):
		# Have to remake the collection from keys because update doesn't work when there are no aligned axes
		# and we can't use aligned axes because all parts of the collection must have aligned axes, which is
		# not necessarily true for the elements of the em collection!!! Obviously we'll want to change this
		# back once those issues are fixed on the NDCube side.
		pairs = []
		for k in self.collection.keys(): pairs.append((k,self.collection.get(k)))
		try:
			pairs.append((modelsequence.meta['ALGORITHM'],modelsequence))
		except KeyError:
			pairs.append((modelsequence.meta['algorithm'],modelsequence))
		self.collection.clear()
		self.collection = NDCollection(pairs)
		try:
			self.collection['models'].append(modelsequence.meta['ALGORITHM'])
		except:
			self.collection['models'].append(modelsequence.meta['algorithm'])

	def compute_dem(self,i,j,logt=None,algorithm=None):
		from scipy.interpolate import interp1d

		if(algorithm is None): algorithm = self.collection['models'][0]
		model = self.collection[algorithm]
		try:
			if(logt is None): logt = model[0].meta['LOGT']
		except KeyError:
			if(logt is None): logt = model[0].meta['logt']
		dem = np.zeros(logt.size)
		for component in model:
			try:
				complogt = component.meta['LOGT']
				compbasis = component.meta['BASIS']
			except KeyError:
				complogt = component.meta['logt']
				compbasis = component.meta['basis']
			dem += component.data[i,j]*interp1d(complogt,compbasis,fill_value=0.0,bounds_error=False)(logt)
		return logt,dem


	def compute_dem_all(self,logt=None,algorithm=None):
		from scipy.interpolate import interp1d

		if(algorithm is None): algorithm = self.collection['models'][0]
		model = self.collection[algorithm]
		if(logt is None): logt = model[0].meta['LOGT']
		dem = np.zeros([model[0].data.shape[0], model[0].data.shape[1], logt.size])
		for component in model:
			complogt = component.meta['LOGT']
			compbasis = component.meta['BASIS']
			dem += np.expand_dims(component.data,-1)*interp1d(complogt,compbasis,fill_value=0.0,bounds_error=False)(logt)
		return logt,dem

	def synthesize_data(self, logts, tresps, algorithm=None, channels=None, ilo=0, ihi=-1, jlo=0, jhi=-1, meta=None):
		from scipy.integrate import trapezoid
		from scipy.interpolate import interp1d
		if(logts[0].size == 1): [logts, tresps] = [[logts],[tresps]]
		ndata = len(logts)
		if(algorithm is None): algorithm = self.collection['models'][0]
		if(channels==None): channels = ['SYNTHDATA'+str(i) for i in range(0,ndata)]
		model = self.collection[algorithm]
		[synthmaps,syntherrs] = [[],[]]
		for i in range(0, ndata):
			synthdata = np.zeros(model[0].data[ilo:ihi,jlo:jhi].shape)
			syntherrs.append(UnknownUncertainty(np.zeros(model[0].data.shape)-1))
			synthmap = copy.deepcopy(model[0])[ilo:ihi,jlo:jhi]
			synthmap.meta['NAXIS'] = 2
			synthmap.meta['ALGORITHM'] = algorithm
			synthmap.meta['CHANNEL'] = channels[i]
			datainterp = interp1d(logts[i], tresps[i], fill_value=0.0, bounds_error=False)
			for component in model:
				try:
					basistemp = component.meta['LOGT']
					basisinterp = interp1d(basistemp, component.meta['BASIS'], fill_value=0.0, bounds_error=False)
				except KeyError:
					basistemp = component.meta['logt']
					basisinterp = interp1d(basistemp, component.meta['basis'], fill_value=0.0, bounds_error=False)
				logt = np.unique(np.hstack([basistemp, logts[i]]))
				coupling = trapezoid(datainterp(logt)*basisinterp(logt),x=logt)
				synthdata += coupling*component.data[ilo:ihi,jlo:jhi]
			synthmap.data[:] = synthdata[:]
			synthmaps.append(synthmap)

		return em_data(synthmaps,syntherrs,logts,tresps,channels=channels)

	def synthesize_map(self, map, logt=None, tresp=None, algorithm=None, channel=None):
		if('CHANNEL' not in map.meta):
			if(channel is None): map.meta['CHANNEL'] = map.meta['DETECTOR']+map.meta['WAVE_STR']
			else: map.meta['CHANNEL'] = channel
		if('SCHEMA' not in map.meta):
			map.meta['LOGT'] = logt
			map.meta['TRESP'] = tresp
			map.meta['SHAPE'] = np.array(map.data.shape)
			map.meta['SCHEMA'] = basic_detector(map.meta)
		if(algorithm is None): algorithm = self.collection['models'][0]
		source = self.collection[algorithm][0].meta['SCHEMA']
		coeffs = np.array([data.data for data in self.collection[algorithm].data]).flatten()
		output_map = copy.deepcopy(map)
		output_map.meta['CHANNEL'] += '_SYNTHETIC'
		output_map.data[::] = (((map.meta['SCHEMA']).fwdop(source))*coeffs).reshape(map.data.shape)
		return output_map


