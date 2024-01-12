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
import ipywidgets as widgets
from IPython.display import display, HTML
import time


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

class dashboard_object(object):
    def __init__(self, em_collection):
        self.emc = em_collection

        [nx,ny] = em_collection.collection[em_collection.collection['models'][0]][0].data.shape
        self.xpt_slider = widgets.IntSlider(min=0, max=nx-1, value=10, step=1, description='xpt', continuous_update=False)
        self.ypt_slider = widgets.IntSlider(min=0, max=ny-1, value=100, step=1, description='ypt', continuous_update=False)
        self.rtemp=widgets.FloatSlider(min=5, max=7, value=5.8, step=0.05, description='rtemp', continuous_update=False)
        self.gtemp=widgets.FloatSlider(min=5, max=7, value=6.1, step=0.05, description='gtemp', continuous_update=False)
        self.btemp=widgets.FloatSlider(min=5, max=7, value=6.4, step=0.05, description='btemp', continuous_update=False)
        self.sigma=widgets.FloatSlider(min=0.025, max=0.5, value=0.125, step=0.01, description='sigma', continuous_update=False)
        self.algorithm=widgets.Dropdown(options=self.emc.collection['models'], description='algorithm', continuous_update=False)

        self.fig = None
        self.ax1 = None
        self.ax2 = None
        self.ax3 = None
        self.ax4 = None
        self.crosshair = None
        self.crosshair_mouseover = None
        self.demplot = None
        self.plotmax = 0.0
        self.count = 0
        self.clicking = False
        self.last_click = 1

        self.red_temp = None
        self.grn_temp = None
        self.blu_temp = None
        self.colorbar = None

        self.last_update_time = 0
        self.click_time = 0

        # Define the custom CSS
        self.custom_css = """
        <style>
        .widget-readout {
            box-shadow: 0px 0px 1px 1px #a9a9a9 inset;
        }
        </style>
        """



    def displays(self):
        ui = widgets.HBox([self.xpt_slider,self.ypt_slider,self.rtemp,self.gtemp,self.btemp,self.sigma,self.algorithm])
        out = widgets.interactive_output(self.widgwrap, {'xpt': self.xpt_slider, 'ypt': self.ypt_slider, 'rtemp': self.rtemp, 'gtemp': self.gtemp, 'btemp': self.btemp, 'sigma': self.sigma, 'algorithm': self.algorithm})
        return ui, out

    def display(self):
        ui, out = self.displays()
        display(ui,out)

    def widgwrap(self, xpt, ypt, rtemp, gtemp, btemp, sigma, algorithm):
        if self.fig is None:
            self.create_figure()
            self.init_figure(xpt, ypt, rtemp, gtemp, btemp, sigma, algorithm)
        else:
            # print(f"widget changed {self.count}, click = {self.clicking}")
            self.count += 1
            if not self.clicking:
                self.update_figure(xpt, ypt, rtemp, gtemp, btemp, sigma, algorithm)


    def create_figure(self):
        # print("Creating Dashboard")
        self.fontsize_prev = plt.rcParams.get('font.size')
        plt.rcParams.update({'font.size':22})

        # Display the custom CSS in the notebook
        HTML(self.custom_css)
        self.fig = plt.figure(constrained_layout=True)
        self.fig.set_size_inches(16, 8)
        spec = self.fig.add_gridspec(ncols=3, nrows=2, width_ratios=[0.1, 0.6, 0.3], height_ratios=[1, 1])

        self.ax1 = self.fig.add_subplot(spec[:, 0])
        self.ax2 = self.fig.add_subplot(spec[:, 1])
        self.ax3 = self.fig.add_subplot(spec[0, 2])
        self.ax4 = self.fig.add_subplot(spec[1, 2])



    def init_figure(self, xpt, ypt, rtemp, gtemp, btemp, sigma, algorithm, gfac=1.0/2.2, plt_emmax=3.0e27):
        # print("Initializing Dashboard")
        [synthchanlogts, synthchantresps] = lognormal_synthetic_channels([rtemp, gtemp, btemp], sigma)
        [cbcoll, cblogts, cbints, cbsigmas] = dem_color_table(synthchanlogts[0])

        [ilo, ihi, jlo, jhi] = [None]*4  # Update as needed
        synthdata = self.emc.synthesize_data(synthchanlogts, synthchantresps, ilo=ilo, ihi=ihi, jlo=jlo, jhi=jhi, algorithm=algorithm)
        self.demimage = np.stack([dat.data for dat in synthdata.data]).T
        colorbar_synthdata = cbcoll.synthesize_data(synthchanlogts, synthchantresps)
        clbimage = np.stack([dat.data for dat in colorbar_synthdata.data]).T

        [i, j] = [xpt, ypt]

        self.demimg = self.ax2.imshow(((np.clip(self.demimage, 0, plt_emmax)/plt_emmax)**gfac).transpose((1, 0, 2)))
        self.crosshair, = self.ax2.plot([i], [j], color='C0', marker='+',           linewidth=8, markersize=20) #, markeredgewidth=3, markeredgecolor='white',)
        self.crosshair_2, = self.ax2.plot([i], [j], color='C1', marker='+',         linewidth=8, markersize=20) #, markeredgewidth=3, markeredgecolor='white',)
        self.crosshair_mouseover, = self.ax2.plot([i], [j], color='C2', marker='+', linewidth=8, markersize=20) #, markeredgewidth=3, markeredgecolor='white',)

        [ptlogt, ptdem] = self.emc.compute_dem(i, j, algorithm=algorithm)
        self.demplot, = self.ax3.plot(10*ptlogt, ptdem/1.0e28,              label=f'Click   at [{i:03}, {j:03}]')
        self.demplot_2, = self.ax3.plot(10*ptlogt, ptdem/1.0e28,            label=f'Click2 at [{i:03}, {j:03}]')
        self.demplot_mouseover, = self.ax3.plot(10*ptlogt, ptdem/1.0e28,    label=f'Mouse at [{i:03}, {j:03}]')
        self.ax3.set_ylim(0, 1.1)
        self.ax3.legend(loc='upper right', fontsize=12)

        plt.suptitle(synthdata[0].meta['algorithm'] + ' inversion at ' + self.emc.data()[0].meta['date-obs'])
        self.ax2.set(title='RGB Composite DEM image')
        self.ax3.set(title='Diff Emission Measure', xlabel='Temperature (dB Kelvin)', ylabel='DEM (Mm \n[$10^9$ cm$^{-3}$]$^2$/dBK)')
        self.ax3.minorticks_on()
        self.ax4.set(title='RGB Composite DEM channel responses', xlabel='Temperature (dB Kelvin)')


        self.red_temp, = self.ax4.plot(10*synthchanlogts[0], synthchantresps[0], 'r')
        self.grn_temp, = self.ax4.plot(10*synthchanlogts[1], synthchantresps[1], 'g')
        self.blu_temp, = self.ax4.plot(10*synthchanlogts[2], synthchantresps[2], 'b')



        self.colorbar = self.ax1.imshow(((clbimage/np.max(clbimage))**gfac), aspect='auto', extent=[cbints[0], cbints[-1], 10*cblogts[0], 10*cblogts[-1]])
        self.ax1.set(title='Color Reference', ylabel='Temperature (dB Kelvin)', xlabel='Channel EM')

        [nx,ny,nz] = self.demimage.shape
        def on_click(event):
            if event.inaxes == self.ax2:
                self.clicking = True
                self.click_time = time.time()
                self.last_click *= -1
                ix, iy = int(event.xdata), int(event.ydata)
                if self.last_click < 0:
                    # Update the xpt and ypt sliders
                    i = self.xpt_slider.value = min(max(ix, 0), nx-1)
                    j = self.ypt_slider.value = min(max(iy, 0), ny-1)
                    self.crosshair.set_data([i], [j])
                    [ptlogt,ptdem] = self.emc.compute_dem(i,j,algorithm=algorithm)
                    self.demplot.set_data(10*ptlogt,ptdem/1.0e28)
                    self.plotmax = max(self.plotmax, np.amax(ptdem/1.0e28))
                    self.demplot.set_label('Click  at '+str([i,j])) #, ylim=(0, 1.1*self.plotmax))
                    self.ax3.legend(loc='upper right', fontsize=12)
                else:
                    # Update the xpt and ypt sliders
                    i = min(max(ix, 0), nx-1)
                    j = min(max(iy, 0), ny-1)
                    self.crosshair_2.set_data([i], [j])
                    [ptlogt,ptdem] = self.emc.compute_dem(i,j,algorithm=algorithm)
                    self.demplot_2.set_data(10*ptlogt,ptdem/1.0e28)
                    self.demplot_2.set_label('Click2 at '+str([i,j])) #, ylim=(0, 1.1*self.plotmax))
                    self.ax3.legend(loc='upper right', fontsize=12)
                self.fig.canvas.draw_idle()
            self.clicking = False
        self.fig.canvas.mpl_connect('button_press_event', on_click)

        def on_mouseover(event):
            current_time = time.time()
            i, j = self.xpt_slider.value, self.ypt_slider.value

            if current_time - self.last_update_time < 1/2:
                return  # Skip the update if less than 1/3 second has passed
            if current_time - self.click_time < 1/4:
                return
            if event.inaxes == self.ax2:
                ix, iy = int(event.xdata), int(event.ydata)

                [ptlogt,ptdem] = self.emc.compute_dem(ix,iy,algorithm=algorithm)
                self.demplot_mouseover.set_data(10*ptlogt,ptdem/1.0e28)

                self.crosshair_mouseover.set_data([ix],[iy])
                self.demplot_mouseover.set_label(f"Mouse at [{ix}, {iy}]")
                self.ax3.legend(loc='upper right', fontsize=12)
                self.fig.canvas.draw_idle()
            else:

                self.crosshair_mouseover.set_data([np.nan],[np.nan])
                self.demplot_mouseover.set_data([np.nan],[np.nan])
                self.demplot_mouseover.set_label("Mouse off chart")
                self.ax3.legend(loc='upper right', fontsize=12)
                self.fig.canvas.draw_idle()

        self.fig.canvas.mpl_connect('motion_notify_event', on_mouseover)



    def update_figure(self, xpt, ypt, rtemp, gtemp, btemp, sigma, algorithm, gfac=1.0/2.2, plt_emmax=3.0e27):
        # Update the plots using the stored handles
        if self.crosshair is not None:
            [i, j] = [xpt, ypt]
            self.crosshair.set_data([i], [j])

        if self.demplot is not None:
            [ptlogt, ptdem] = self.emc.compute_dem(i, j, algorithm=algorithm)
            self.demplot.set_data(10*ptlogt, ptdem/1.0e28)

            self.plotmax = max(self.plotmax, np.amax(ptdem/1.0e28))

        [synthchanlogts, synthchantresps] = lognormal_synthetic_channels([rtemp, gtemp, btemp], sigma)
        if self.red_temp is not None:
            self.red_temp.set_data(10*synthchanlogts[0], synthchantresps[0])
            self.grn_temp.set_data(10*synthchanlogts[1], synthchantresps[1])
            self.blu_temp.set_data(10*synthchanlogts[2], synthchantresps[2])

        if self.colorbar is not None:
            [cbcoll, cblogts, cbints, cbsigmas] = dem_color_table(synthchanlogts[0])
            colorbar_synthdata = cbcoll.synthesize_data(synthchanlogts, synthchantresps)
            clbimage = np.stack([dat.data for dat in colorbar_synthdata.data]).T
            self.colorbar.set_data(((clbimage/np.max(clbimage))**gfac))
            self.colorbar.set_extent([cbints[0], cbints[-1], 10*cblogts[0], 10*cblogts[-1]])

        if self.demimg is not None:
            [ilo, ihi, jlo, jhi] = [None]*4  # Update as needed
            synthdata = self.emc.synthesize_data(synthchanlogts, synthchantresps, ilo=ilo, ihi=ihi, jlo=jlo, jhi=jhi, algorithm=algorithm)
            demimage = np.stack([dat.data for dat in synthdata.data]).T
            self.demimg.set_data(((np.clip(demimage, 0, plt_emmax)/plt_emmax)**gfac).transpose((1, 0, 2)))

        self.fig.canvas.draw()

