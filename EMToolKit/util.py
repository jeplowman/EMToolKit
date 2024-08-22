import copy
import numpy as np
from sunpy.map import Map
from ndcube import NDCube, NDCubeSequence, NDCollection
from astropy.coordinates import SkyCoord
from astropy.nddata import StdDevUncertainty
import matplotlib.pyplot as plt

def lognormal_synthetic_channels(temps,sigmas=0.1,logt=None,nt=81):

    if(np.size(sigmas) == 1): sigmas = sigmas+0.0*np.array(temps)
    if(logt is None): logt = np.linspace(min(temps)-2.5*np.max(sigmas),max(temps)+2.5*np.max(sigmas),nt)

    tresps = []
    for i in range(0,len(temps)):
        tresps.append(np.exp(-0.5*(logt-temps[i])**2.0/(sigmas[i])**2.0))
    logts = [logt for temp in temps]

    return logts,tresps

def triangle_basis(logt):
    nt = len(logt)
    basislogt = np.linspace(np.min(logt),np.max(logt),2*(nt-1)+1)
    [basislogts,bases] = [[],[]]
    for i in range(0,nt):
        basis = (basislogt == logt[i]).astype(np.float64)
        if(i > 0): basis += (basislogt-logt[i-1])*(basislogt < logt[i])*(basislogt > logt[i-1])/(logt[i]-logt[i-1])
        if(i < nt-1): basis += (logt[i+1]-basislogt)*(basislogt < logt[i+1])*(basislogt > logt[i])/(logt[i+1]-logt[i])
        bases.append(basis)
        basislogts.append(basislogt)
    return basislogts,bases

def calc_resids(synthdata, em_collection):# Calculate the residuals and Chi squared:
    ndata = len(synthdata)
    resids = []
    datasequence = em_collection.data()
    chi2 = 0
    [nx,ny] = datasequence[0].data.shape
    for seq in datasequence: [nx,ny] = [np.min([seq.data.shape[0],nx]),np.min([seq.data.shape[1],ny])]
    for i in range(0,ndata):
        exptime = datasequence[i].meta['exptime']
        nx = np.min([synthdata[i].data.shape[0],datasequence[i].data.shape[0],datasequence[i].uncertainty.array.shape[0]])
        ny = np.min([synthdata[i].data.shape[1],datasequence[i].data.shape[1],datasequence[i].uncertainty.array.shape[1]])
        resids.append(((exptime*synthdata[i].data[0:nx,0:ny]-datasequence[i].data[0:nx,0:ny])/datasequence[i].uncertainty.array[0:nx,0:ny])**2)
        chi2 += np.mean(resids)/ndata
    return resids, chi2

def plot_resids(em_collection,resids,algorithm, figsize=None):
    fig = plt.figure(figsize=figsize)
    plt.suptitle('Residuals for '+algorithm)
    for i in range(0,6):
        ax1 = fig.add_subplot(2,3,i+1)
        ax1.imshow(resids[i],vmin=0,vmax=5)
        ax1.set(title=em_collection.data()[i].meta['channel'])


import os.path
def list_fits_files(sdo_data_dir, find=None):
    paths = [os.path.join(sdo_data_dir, path) for path in os.listdir(sdo_data_dir)
            if os.path.isfile(os.path.join(sdo_data_dir, path)) and ".fits" in path]
    if find is not None:
        paths = [path for path in paths if find in path or find.casefold() in path.casefold()]
    return paths