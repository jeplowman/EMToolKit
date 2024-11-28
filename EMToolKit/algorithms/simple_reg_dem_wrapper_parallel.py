import numpy as np
from sunpy.map import Map
from ndcube import NDCube, NDCubeSequence, NDCollection
from astropy.coordinates import SkyCoord
from astropy.nddata import StdDevUncertainty
from EMToolKit.algorithms.simple_reg_dem_parallel import simple_reg_dem
import os.path
import pickle
import time
import EMToolKit.EMToolKit as emtk


# Need to implement passing the wrapargs to the init routines...
def simple_reg_dem_wrapper(datasequence,wrapargs=None):
    nc = len(datasequence)
    prepend = wrapargs.get('prepend',"")
    [nx,ny] = datasequence[0].data.shape
    for seq in datasequence: [nx,ny] = [np.min([seq.data.shape[0],nx]),np.min([seq.data.shape[1],ny])]
    logt = datasequence[0].meta['logt']
    datacube = np.zeros([nx,ny,nc])
    errscube = np.zeros([nx,ny,nc])
    tresps = np.zeros([logt.size,nc])
    exptimes = np.zeros(nc)
    for i in range(0,nc):
        datacube[:,:,i] = datasequence[i].data[0:nx,0:ny]
        errscube[:,:,i] = datasequence[i].uncertainty.array[0:nx,0:ny]
        tresps[:,i] = datasequence[i].meta['tresp']
        exptimes[i] = datasequence[i].meta['exptime']

    coeffs,chi2 = simple_reg_dem(datacube,errscube,exptimes,logt,tresps)
    # Simple_reg_dem puts the temperature axis last. Transpose so it's the first:
    coeffs = coeffs.transpose(np.roll(np.arange(coeffs.ndim),1))

    nt = logt.size
    wcs = datasequence[0].wcs
    basislogt = np.linspace(np.min(logt),np.max(logt),2*(nt-1)+1)
    [logts,bases] = [[],[]]
    for i in range(0,nt):
        basis = (basislogt == logt[i]).astype(np.float64)
        if(i > 0): basis += (basislogt-logt[i-1])*(basislogt < logt[i])*(basislogt > logt[i-1])/(logt[i]-logt[i-1])
        if(i < nt-1): basis += (logt[i+1]-basislogt)*(basislogt < logt[i+1])*(basislogt > logt[i])/(logt[i+1]-logt[i])
        bases.append(basis)
        logts.append(basislogt)
    return list(coeffs),logts,bases,wcs,prepend+'simple_reg_dem'



def autoloading_simple_reg_dem_wrapper(datasequence, data_dir="../", recalc_simple=False, wrapargs={}):
    pk_file = os.path.join(data_dir, 'simple_reg_demsequence.pkl')

    if os.path.exists(pk_file) and not recalc_simple:
        print('Loading simple_reg_demsequence from', pk_file)
        with open(pk_file, 'rb') as file:
            (simple_reg_demsequence, simpl_out) = pickle.load(file)
    else:
        print("Calculating simple_reg_demsequence from scratch in parallel...", end="")
        tstart=time.time()
        simpl_out = simple_reg_dem_wrapper(datasequence, wrapargs)
        print('Done! Simple method took',time.time()-tstart)
        simple_reg_demsequence = emtk.dem_model(*simpl_out, simple_reg_dem_wrapper)
        if not os.path.exists(data_dir):
            os.makedirs(data_dir)
        with open(pk_file, 'wb') as file:
            pickle.dump((simple_reg_demsequence, simpl_out), file)
    return simple_reg_demsequence, simpl_out
