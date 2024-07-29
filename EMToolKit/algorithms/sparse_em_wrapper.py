import copy
import numpy as np
from sunpy.map import Map
from ndcube import NDCube, NDCubeSequence, NDCollection
from astropy.coordinates import SkyCoord
from astropy.nddata import StdDevUncertainty
from EMToolKit.algorithms.sparse_em import sparse_em_init, sparse_em_solve
import os.path
import pickle
import time
import EMToolKit.EMToolKit as emtk

# Need to implement passing wrapargs to the sparse_em routines...
def sparse_em_wrapper(datasequence, wrapargs=None):
	nc = len(datasequence)
	[nx,ny] = datasequence[0].data.shape
	for seq in datasequence:
		[nx,ny] = [np.min([seq.data.shape[0],nx]),np.min([seq.data.shape[1],ny])]
	logt = datasequence[0].meta['logt']
	datacube = np.zeros([nx,ny,nc])
	errscube = np.zeros([nx,ny,nc])
	[trlogts,tresps] = [[],[]]
	exptimes = np.zeros(nc)
	for i in range(0,nc):
		datacube[:,:,i] = datasequence[i].data[0:nx,0:ny]
		errscube[:,:,i] = datasequence[i].uncertainty.array[0:nx,0:ny]
		tresps.append(datasequence[i].meta['tresp'])
		trlogts.append(datasequence[i].meta['logt'])
		exptimes[i] = datasequence[i].meta['exptime']

	[Dict, lgtaxis, basis_funcs, bases_sigmas] = sparse_em_init(trlogts, tresps, differential=True)
	coeffs, zmax, status = sparse_em_solve(datacube, errscube, exptimes, Dict)
	coeffs = coeffs.transpose(np.roll(np.arange(coeffs.ndim),1))

	[nchannels,nb] = Dict.shape
	wcs = datasequence[0].wcs
	logts = nb*[lgtaxis]
	return list(coeffs),logts,list(basis_funcs),wcs,'sparse_em'

def autoloading_sparse_em_wrapper(datasequence, data_dir="../.data/default", recalc_sparse=False):
	if not os.path.exists(data_dir):
		os.makedirs(data_dir)

	pk_file = os.path.join(data_dir, 'sparse_em_demsequence.pkl')

	if os.path.exists(pk_file) and not recalc_sparse:
		with open(pk_file, 'rb') as file:
			(sparse_em_demsequence, spars_out) = pickle.load(file)
	else:
		tstart=time.time()
		spars_out = sparse_em_wrapper(datasequence)
		sparse_em_demsequence = emtk.dem_model(*spars_out, sparse_em_wrapper)
		print('Sparse method took',time.time()-tstart)
		with open(pk_file, 'wb') as file:
			pickle.dump((sparse_em_demsequence, spars_out), file)
	return sparse_em_demsequence, spars_out