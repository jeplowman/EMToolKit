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
from scipy.interpolate import make_interp_spline
import scipy.special
from matplotlib import rcParams


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

# class dashboard_object(object):
#     def __init__(self,em_collection):
#         self.emc = em_collection

#     def widgwrap(self, xpt, ypt, rtemp, gtemp, btemp, sigma, algorithm):
#         dashboard_figure(self.emc, plotpoint=[xpt,ypt], temperatures=[rtemp,gtemp,btemp], sigmas=sigma, algorithm=algorithm)

class dashboard_object(object):
    def __init__(self, em_collection):
        self.emc = em_collection

        [nx,ny] = em_collection.collection[em_collection.collection['models'][0]][0].data.shape
        # self.xpt_slider = widgets.IntSlider(min=0, max=nx-1, value=10, step=1, description='xpt', continuous_update=False)
        # self.ypt_slider = widgets.IntSlider(min=0, max=ny-1, value=100, step=1, description='ypt', continuous_update=False)
        self.rtemp=widgets.FloatSlider(min=5, max=7, value=5.8, step=0.05, description='rtemp', continuous_update=False)
        self.gtemp=widgets.FloatSlider(min=5, max=7, value=6.1, step=0.05, description='gtemp', continuous_update=False)
        self.btemp=widgets.FloatSlider(min=5, max=7, value=6.4, step=0.05, description='btemp', continuous_update=False)
        self.sigma=widgets.FloatSlider(min=0.025, max=0.5, value=0.125, step=0.01, description='sigma', continuous_update=False)
        self.algorithm=widgets.Dropdown(options=self.emc.collection['models'], description='algorithm', continuous_update=False)
        self.btn_draw_curve = widgets.Button(description="Draw Curve")
        self.btn_finish_lines = widgets.Button(description="Reset Lines")

        self.curve_points = []
        self.drawing = False




        self.fig = None
        self.ax1 = None
        self.ax2 = None
        self.ax3 = None
        self.ax4 = None
        self.crosshair = None
        self.crosshair_mouseover = None
        self.bezline = None
        self.demplot = None
        self.plotmax = 0.0
        self.count = 0
        self.clicking = False
        self.last_click = 1

        self.red_temp = None
        self.grn_temp = None
        self.blu_temp = None
        self.colorbar = None
        self.the_algorithm = None
        self.crosshairs = []
        self.crosshairs_bezier = []
        self.demlines = []
        self.legend = None

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
        ui0 = widgets.HBox([self.rtemp,self.gtemp,self.btemp,self.sigma,self.algorithm])
        ui1 = widgets.HBox([self.btn_draw_curve, self.btn_finish_lines])
        ui = widgets.VBox([ui0,ui1])
        out = widgets.interactive_output(self.widgwrap, {'rtemp': self.rtemp, 'gtemp': self.gtemp, 'btemp': self.btemp, 'sigma': self.sigma, 'algorithm': self.algorithm})
        return ui, out

    def display(self):
        ui, out = self.displays()
        display(ui,out)

    def widgwrap(self, rtemp, gtemp, btemp, sigma, algorithm):
        if self.fig is None:
            self.create_figure()
            self.init_figure( rtemp, gtemp, btemp, sigma, algorithm)
        else:
            self.count += 1
            self.update_figure( rtemp, gtemp, btemp, sigma, algorithm)


    def create_figure(self):
        # print("Creating Dashboard")
        self.fontsize_prev = plt.rcParams.get('font.size')
        plt.rcParams.update({'font.size':18})

        # Display the custom CSS in the notebook
        HTML(self.custom_css)
        self.fig = plt.figure(constrained_layout=True)
        self.fig.set_size_inches(16, 8)
        spec = self.fig.add_gridspec(ncols=3, nrows=2, width_ratios=[0.1, 0.6, 0.3], height_ratios=[1, 1])

        self.ax1 = self.fig.add_subplot(spec[:, 0])
        self.ax2 = self.fig.add_subplot(spec[:, 1])
        self.ax3 = self.fig.add_subplot(spec[0, 2])
        self.ax4 = self.fig.add_subplot(spec[1, 2])

    def update_legend(self):
        # Hide the legend based on the condition
        self.legend = self.ax3.legend(loc='upper right', fontsize=12,bbox_to_anchor=(1, 1))

        # if self.count > 5:
        #     if self.legend is not None:
        #         self.legend.set_visible(False)
        # else:
        #     if self.legend is not None:
        #         self.legend.set_visible(True)

        self.fig.canvas.draw_idle()

    def init_dem_line(self, ix, iy):
        NC = self.count
        self.crosshairs.append(self.ax2.plot([ix], [iy], marker='+', color=f"C{NC}", markersize=25)[0])
        [ptlogt, ptdem] = self.emc.compute_dem(ix, iy, algorithm=self.the_algorithm)
        self.count += 1
        thelabel = f'Click {self.count} at [{ix:03}, {iy:03}]' if self.count < 6 else None
        self.demlines.append(self.ax3.plot(10*ptlogt, ptdem/1.0e28, color=f"C{NC}", label=thelabel)[0])
        self.update_legend()

    def init_mouseover_line(self):
        NC = self.count
        self.crosshair_mouseover, = self.ax2.plot([], [], color='purple', marker='+', markersize=25)
        # [ptlogt, ptdem] = self.emc.compute_dem(0, 0, algorithm=self.the_algorithm)
        self.demplot_mouseover, = self.ax3.plot([],[], color='purple', ls="--", zorder=10000,  label=f"Mouse off chart")
        self.update_legend()





    def init_figure(self, rtemp, gtemp, btemp, sigma, algorithm, gfac=1.0/2.2, plt_emmax=3.0e27):
        # print("Initializing Dashboard")
        self.the_algorithm = algorithm

        [synthchanlogts, synthchantresps] = lognormal_synthetic_channels([rtemp, gtemp, btemp], sigma)
        [cbcoll, cblogts, cbints, cbsigmas] = dem_color_table(synthchanlogts[0])

        [ilo, ihi, jlo, jhi] = [None]*4  # Update as needed
        synthdata = self.emc.synthesize_data(synthchanlogts, synthchantresps, ilo=ilo, ihi=ihi, jlo=jlo, jhi=jhi, algorithm=self.the_algorithm)
        self.demimage = np.stack([dat.data for dat in synthdata.data]).T
        colorbar_synthdata = cbcoll.synthesize_data(synthchanlogts, synthchantresps)
        clbimage = np.stack([dat.data for dat in colorbar_synthdata.data]).T
        self.demimg = self.ax2.imshow(((np.clip(self.demimage, 0, plt_emmax)/plt_emmax)**gfac).transpose((1, 0, 2)))

        [ptlogt, ptdem] = self.emc.compute_dem(0,0, algorithm=self.the_algorithm)
        self.init_mouseover_line()
        self.ax3.set_ylim(0, 1.1)
        self.ax3.set_xlim(np.min(10*ptlogt), np.max(10*ptlogt))
        self.update_legend()

        try:
            plt.suptitle(synthdata[0].meta['algorithm'] + ' inversion at ' + self.emc.data()[0].meta['date-obs'])
        except KeyError:
            plt.suptitle(synthdata[0].meta['ALGORITHM'] + ' inversion at ' + self.emc.data()[0].meta['date-obs'])

        self.ax2.set(title='RGB Composite DEM image')
        self.ax3.set(title='Diff Emission Measure', xlabel='Temperature (dB Kelvin)', ylabel='DEM (Mm \n[$10^9$ cm$^{-3}$]$^2$/dBK)')
        self.ax3.minorticks_on()
        self.ax4.set(title='RGB Composite DEM channel responses', xlabel='Temperature (dB Kelvin)')


        self.red_temp, = self.ax4.plot(10*synthchanlogts[0], synthchantresps[0], 'r')
        self.grn_temp, = self.ax4.plot(10*synthchanlogts[1], synthchantresps[1], 'g')
        self.blu_temp, = self.ax4.plot(10*synthchanlogts[2], synthchantresps[2], 'b')

        self.colorbar = self.ax1.imshow(((clbimage/np.max(clbimage))**gfac), aspect='auto', extent=[cbints[0], cbints[-1], 10*cblogts[0], 10*cblogts[-1]])
        self.ax1.set(title='Color Reference', ylabel='Temperature (dB Kelvin)', xlabel='Channel EM')

        self.init_interactivity()

    def init_interactivity(self):

        [nx,ny,nz] = self.demimage.shape
        def on_click(event):
            if event.inaxes == self.ax2:
                ix, iy = int(event.xdata), int(event.ydata)
                i = min(max(ix, 0), nx-1)
                j = min(max(iy, 0), ny-1)

                if self.drawing:
                    self.update_bezier_curve(i, j)  # Function to draw/update the Bezier curve
                # else:
                self.init_dem_line(i, j)

            self.update_legend()

        self.fig.canvas.mpl_connect('button_press_event', on_click)

        def on_mouseover(event):

            if self.demplot_mouseover is None:
                self.init_mouseover_line()

            if event.inaxes == self.ax2:
                ix, iy = int(event.xdata), int(event.ydata)

                [ptlogt,ptdem] = self.emc.compute_dem(ix,iy,algorithm=self.the_algorithm)
                self.demplot_mouseover.set_data(10*ptlogt,ptdem/1.0e28)
                self.crosshair_mouseover.set_data([ix],[iy])
                self.demplot_mouseover.set_label(f"Mouse at [{ix}, {iy}]")
                self.update_legend()
            else:
                self.crosshair_mouseover.set_data([np.nan],[np.nan])
                self.demplot_mouseover.set_data([np.nan],[np.nan])
                self.demplot_mouseover.set_label("Mouse off chart")
                self.update_legend()

        self.fig.canvas.mpl_connect('motion_notify_event', on_mouseover)

        def on_draw_curve_clicked(b):
            print("Drawing curve")
            if b.description == "Draw Curve":
                self.drawing = True
                self.curve_points = []  # Reset the points
                b.description = "Stop Drawing"
                b.button_style = 'success'  # Red color
            else:
                self.drawing = False
                b.description = "Draw Curve"
                b.button_style = ''  # Default color  # Default color

        def on_reset_lines_clicked(b):
            for line in self.ax3.lines:
                line.remove()
            for crosshair in self.crosshairs:
                if crosshair is not None and crosshair.axes is not None:
                    crosshair.remove()
            if self.bezline is not None:
                self.bezline.set_data([np.nan], [np.nan])

            self.crosshairs = []

            self.count = 0
            self.ax3.set_prop_cycle(rcParams['axes.prop_cycle'])

            self.init_mouseover_line()

            self.update_legend()

            print("Reset Lines")



        self.btn_draw_curve.on_click(on_draw_curve_clicked)
        self.btn_finish_lines.on_click(on_reset_lines_clicked)


    def update_bezier_curve(self, i, j):
        self.curve_points.append((i, j))
        self.crosshairs.append(self.ax2.plot([i], [j], marker='o', markeredgecolor=f"cyan", markersize=10)[0])
        if len(self.curve_points) < 2:
            return


        def compute_bezier_points(control_points, num_points=100):
            n = len(control_points) - 1
            t = np.linspace(0, 1, num_points)
            curve_points = np.zeros((num_points, 2))
            for i in range(n + 1):
                binomial_coeff = scipy.special.comb(n, i)
                # We need to ensure that the shapes are compatible for broadcasting
                # The np.newaxis helps in aligning the shapes for the operation
                curve_points += binomial_coeff * (t**i)[:, np.newaxis] * ((1 - t)**(n - i))[:, np.newaxis] * np.array(control_points[i])
            return curve_points

        # Compute the Bezier curve points
        bezier_points = compute_bezier_points(self.curve_points, num_points=100)

        # Clear the previous curve and draw a new one
        if self.bezline is not None:
            self.bezline.set_data(bezier_points[:, 0], bezier_points[:, 1])
        else:
            self.bezline, = self.ax2.plot(bezier_points[:, 0], bezier_points[:, 1], 'r-')
        self.fig.canvas.draw_idle()



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

        self.fig.canvas.draw_idle()

