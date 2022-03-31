import copy
import numpy as np
from sunpy.map import Map
from ndcube import NDCube, NDCubeSequence, NDCollection
from astropy.coordinates import SkyCoord
from astropy.nddata import StdDevUncertainty
import EMToolKit.EMToolKit as emtk
import matplotlib.pyplot as plt
from EMToolKit.util import lognormal_synthetic_channels, triangle_basis

def dem_color_table(ctlogts,sigmin=0.1,sigmax=0.1,intmin=0.0,intmax=1.0,n=81):

    from astropy import wcs
    [basislogts,bases] = triangle_basis(ctlogts)
    ctints = np.linspace(intmin,intmax,n)
    ctsigmas = np.linspace(sigmin,sigmax,n)
    logt = ctlogts

    nt = len(logt)
    nctlt = len(ctlogts)
    ns = len(ctsigmas)
    clrtab = np.zeros([nt,ns,nctlt])
    for i in range(0,ns):
        for j in range(0,nctlt):
            clrtab[:,i,j] = ctints[i]*np.exp(-0.5*(logt-ctlogts[j])**2.0/(ctsigmas[i])**2.0)
        
    ctmodel = emtk.dem_model(clrtab,basislogts,bases,wcs.WCS(naxis=2),'Color Table',None,meta={},wrapargs=None)
    ctcoll = emtk.em_collection(None)
    ctcoll.add_model(ctmodel)
    return ctcoll,ctlogts,ctints,ctsigmas


def dashboard_figure(em_collection,plotpoint=None,temperatures=[5.8,6.1,6.4], gfac=1.0/2.2,
					sigmas=0.1, cropto = [None]*4, plt_emmax=3.0e27, algorithm=None):
    [synthchanlogts,synthchantresps] = lognormal_synthetic_channels(temperatures,sigmas)
    [cbcoll,cblogts,cbints,cbsigmas] = dem_color_table(synthchanlogts[0])

    [ilo,ihi,jlo,jhi] = cropto
    synthdata = em_collection.synthesize_data(synthchanlogts,synthchantresps,
    					ilo=ilo,ihi=ihi,jlo=jlo,jhi=jhi,algorithm=algorithm)
    demimage = np.stack([dat.data for dat in synthdata.data]).T

    colorbar_synthdata = cbcoll.synthesize_data(synthchanlogts,synthchantresps)
    clbimage = np.stack([dat.data for dat in colorbar_synthdata.data]).T

    if(plotpoint is None): plotpoint = np.round(0.5*np.array(demimage.shape)).astype(np.int32)
    [i,j] = plotpoint[0:2]

    fontsize_prev = plt.rcParams.get('font.size')
    plt.rcParams.update({'font.size':24})

    fig = plt.figure(constrained_layout=True)
    plt.suptitle(synthdata[0].meta['algorithm'] + ' inversion at ' + em_collection.data()[0].meta['date-obs'])
    spec = fig.add_gridspec(ncols = 3, nrows=2,width_ratios = [0.1,0.6,0.3],height_ratios=[1,1])

    ax1 = fig.add_subplot(spec[:,0])
    ax2 = fig.add_subplot(spec[:,1])
    ax3 = fig.add_subplot(spec[0,2])
    ax4 = fig.add_subplot(spec[1,2])

    [ptlogt,ptdem] = em_collection.compute_dem(i,j,algorithm=algorithm)

    ax2.imshow(((np.clip(demimage,0,plt_emmax)/plt_emmax)**gfac).transpose((1,0,2)))
    ax2.plot([i],[j],color='white',marker='+',linewidth=4,markersize=20)
    ax2.set(title='RGB Composite DEM image')
    ax3.plot(10*ptlogt,ptdem/1.0e28)
    ax3.set(title='DEM at '+str([i,j]),xlabel='Temperature (dB Kelvin)',ylabel='DEM (Mm [$10^9$ cm$^{-3}$]$^2$/dBK)')
    ax3.minorticks_on()
    ax4.plot(10*synthchanlogts[0],synthchantresps[0],'r')
    ax4.plot(10*synthchanlogts[1],synthchantresps[1],'g')
    ax4.plot(10*synthchanlogts[2],synthchantresps[2],'b')
    ax4.set(title='RGB Composite DEM channel responses',xlabel='Temperature (dB Kelvin)')
    ax1.imshow(((clbimage/np.max(clbimage))**gfac),aspect='auto',extent=[cbints[0],cbints[-1],10*cblogts[0],10*cblogts[-1]])
    ax1.set(title='Color Reference',ylabel='Temperature (dB Kelvin)',xlabel='Channel EM')

    plt.rcParams.update({'font.size':fontsize_prev})