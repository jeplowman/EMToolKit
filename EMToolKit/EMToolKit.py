import copy
import numpy as np
from sunpy.map import Map
from ndcube import NDCube, NDCubeSequence, NDCollection
from astropy.coordinates import SkyCoord
from astropy.nddata import StdDevUncertainty

def em_data(maps,errs,logts,tresps,channels=None):
    cubes = []
    for i in range(0,len(maps)):
        mapi = copy.copy(maps[i])
        if('channel' not in mapi.meta):
            if(channels is None): mapi.meta['channel'] = mapi.meta['detector']+mapi.meta['wave_str']
            else: mapi.meta['channel'] = channels[i]
        mapi.meta['logt'] = logts[i]
        mapi.meta['tresp'] = tresps[i]
        cubes.append(NDCube(mapi,uncertainty=errs[i]))
    return NDCubeSequence(cubes)
    
def dem_model(coeffs,logts,bases,wcs,algorithm,wrapper,meta={},wrapargs=None):

    nd = len(coeffs)
    dem_sequence = []
    for i in range(0,nd):
        logt0 = np.median(logts[i][np.where(bases[i]==np.max(bases[i]))])
        metai = {'basis':bases[i],'logt':logts[i],'logt0':logt0,'algorithm':algorithm,'wrapper':wrapper,'wrapargs':wrapargs}
        dem_sequence.append(NDCube(coeffs[i],wcs=wcs,meta={**meta,**metai}))
    return NDCubeSequence(dem_sequence,meta={'algorithm':algorithm,'wrapper':wrapper,'wrapargs':wrapargs})
    
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
		pairs.append((modelsequence.meta['algorithm'],modelsequence))
		self.collection.clear()
		self.collection = NDCollection(pairs)
		self.collection['models'].append(modelsequence.meta['algorithm'])
		