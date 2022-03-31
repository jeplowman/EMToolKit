import copy
import numpy as np
from sunpy.map import Map
from ndcube import NDCube, NDCubeSequence, NDCollection
from astropy.coordinates import SkyCoord
from astropy.nddata import StdDevUncertainty

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