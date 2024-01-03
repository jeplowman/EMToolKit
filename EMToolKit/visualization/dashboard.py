from __future__ import print_function
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
        
    ctmodel = emtk.dem_model(clrtab,basislogts,bases,wcs.WCS(naxis=2),'Color Table',None,meta=dummy_meta(nt,ns),wrapargs=None)
    ctcoll = emtk.em_collection(None)
    ctcoll.add_model(ctmodel)
    return ctcoll,ctlogts,ctints,ctsigmas

def dummy_meta(nx,ny):
	return {'CDELT1':1,'CDELT2':1,'CROTA':0,'CRPIX1':1,'CRPIX2':1,'CRVAL1':0,'CRVAL2':0,'NAXIS1':nx,'NAXIS2':ny,
			'CUNIT1':'arcsec','CUNIT2':'arcsec','CTYPE1':'HPLN-TAN','CTYPE2':'HPLT-TAN','DATE-OBS':'1980-01-01T00:00:00.000'}

# May need installation, for example some of the following
# pip install ipywidgets
# conda install -c conda-forge ipywidgets
# conda install -n base -c conda-forge widgetsnbextension
# conda install -n py36 -c conda-forge ipywidgets
# see https://ipywidgets.readthedocs.io/en/latest/user_install.html

from ipywidgets import interact, interactive, fixed, interact_manual
import ipywidgets as widgets

class dashboard_object(object):
    def __init__(self,em_collection):
        self.emc = em_collection
        
    def widgwrap(self, xpt, ypt, rtemp, gtemp, btemp, sigma, algorithm):
        dashboard_figure(self.emc, plotpoint=[xpt,ypt], temperatures=[rtemp,gtemp,btemp], sigmas=sigma, algorithm=algorithm)

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
	plt.suptitle(synthdata[0].meta['ALGORITHM'] + ' inversion at ' + em_collection.data()[0].meta['DATE-OBS'])
	spec = fig.add_gridspec(ncols = 3, nrows=2,width_ratios = [0.1,0.6,0.3],height_ratios=[1,1])

	ax1 = fig.add_subplot(spec[:,0])
	ax2 = fig.add_subplot(spec[:,1])
	ax3 = fig.add_subplot(spec[0,2])
	ax4 = fig.add_subplot(spec[1,2])

	[ptlogt,ptdem] = em_collection.compute_dem(i,j,algorithm=algorithm)

	ax2.imshow(((np.clip(demimage,0,plt_emmax)/plt_emmax)**gfac).transpose((1,0,2)))
	ax2.plot([i],[j],color='white',marker='+',linewidth=4,markersize=20)
	ax2.set(title='RGB Composite DEM image')
	ax3.semilogx(10**ptlogt,ptdem/1.0e28)
	ax3.set(title='DEM at '+str([i,j]),xlabel='Temperature (Kelvin)',ylabel='DEM (Mm [$10^9$ cm$^{-3}$]$^2$/dBK)')
	ax3.minorticks_on()
	ax4.semilogx(10**synthchanlogts[0],synthchantresps[0],'r')
	ax4.semilogx(10**synthchanlogts[1],synthchantresps[1],'g')
	ax4.semilogx(10**synthchanlogts[2],synthchantresps[2],'b')
	ax4.set(title='RGB Composite DEM channel responses',xlabel='Temperature (Kelvin)')
	ax1.imshow(((clbimage/np.max(clbimage))**gfac),aspect='auto',extent=[cbints[0],cbints[-1],1*cblogts[0],1*cblogts[-1]])

	fig.canvas.draw()

	labels = [str(round(10**float(item.get_text())/1.0e4)/100) for item in ax1.get_yticklabels()]

	ax1.set_yticklabels(labels)
	ax1.set(title='Color Reference',ylabel='Temperature (MKelvin)',xlabel='Channel EM')

	plt.show()

	plt.rcParams.update({'font.size':fontsize_prev})
