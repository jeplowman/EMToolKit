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
from IPython.display import display, HTML, clear_output
import time
from scipy.interpolate import make_interp_spline
import scipy.special
from matplotlib import rcParams
from scipy.interpolate import splprep, splev, UnivariateSpline

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
	return {'CDELT1':1,'CDELT2':1,'CROTA':0,'CRPIX1':1,'CRPIX2':1,'CRVAL1':0,'CRVAL2':0,'NAXIS1':nx,'NAXIS2':ny,'HGLT_OBS':0.5*np.pi,'HGLN_OBS':0.0,'DSUN_OBS':1.49e11,
			'CUNIT1':'arcsec','CUNIT2':'arcsec','CTYPE1':'HPLN-TAN','CTYPE2':'HPLT-TAN','DATE-OBS':'1980-01-01T00:00:00.000','RSUN_REF':696000000}

# May need installation, for example some of the following
# pip install ipywidgets
# conda install -c conda-forge ipywidgets
# conda install -n base -c conda-forge widgetsnbextension
# conda install -n py36 -c conda-forge ipywidgets
# see https://ipywidgets.readthedocs.io/en/latest/user_install.html

from ipywidgets import interact, interactive, fixed, interact_manual
import ipywidgets as widgets


class dashboard_object(object):



    def __init__(self, em_collection, **kwargs):
        self.emc = em_collection
        self.first = em_collection.collection[em_collection.collection['models'][0]][0]
        
        rt0,gt0,bt0 = kwargs.get('rtemp',5.6),kwargs.get('gtemp',6.1),kwargs.get('btemp',6.6)
        sg0 = 0.5*np.mean(np.sort([rt0,gt0,bt0])[1:]-np.sort([rt0,gt0,bt0])[0:-1])

        [nx,ny] = em_collection.collection[em_collection.collection['models'][0]][0].data.shape
        # self.xpt_slider = widgets.IntSlider(min=0, max=nx-1, value=10, step=1, description='xpt', continuous_update=False)
        # self.ypt_slider = widgets.IntSlider(min=0, max=ny-1, value=100, step=1, description='ypt', continuous_update=False)
        self.rtemp=widgets.FloatSlider(min=5, max=7, value=rt0, step=0.05, description='rtemp', continuous_update=False)
        self.gtemp=widgets.FloatSlider(min=5, max=7, value=gt0, step=0.05, description='gtemp', continuous_update=False)
        self.btemp=widgets.FloatSlider(min=5, max=7, value=bt0, step=0.05, description='btemp', continuous_update=False)
        self.sigma=widgets.FloatSlider(min=0.025, max=0.5, value=sg0, step=0.01, description='sigma', continuous_update=False)
        self.rng=widgets.FloatRangeSlider(min=55, max=75, value=(58, 68), step=0.5, description='PlotRange', continuous_update=False)
        self.algorithm=widgets.Dropdown(options=self.emc.collection['models'], description='algorithm', continuous_update=False)
        self.normalization=widgets.Dropdown(options=['max', 'area', 'none'], description='norm', continuous_update=True)
        self.init_buttons()
        self.xsize,self.ysize=kwargs.get('xsize',14),kwargs.get('ysize',10)
        self.fontsize = kwargs.get('fontsize',18)

        # self.slice_type = "bezier"
        self.slice_type=widgets.Dropdown(options=["spline","bezier"], description='slice type', continuous_update=False)
        self.mouseover = widgets.Checkbox( value=True, description='mouseover')
        self.tick_spacing = widgets.IntSlider(min=5, max=100, value=50, step=5, description='spacing', continuous_update=False)
        self.tick_spacing_value = 50
        self.control_points = []
        self.drawing = True
        self.the_normalization = "none"
        self.demplot_mouseover_vert = None
        self.demplot_mouseover = None
        self.fig = None
        self.ax1 = None
        self.ax2 = None
        self.ax3 = None
        self.ax4 = None
        self.crosshair = None
        self.crosshair_mouseover = None
        self.slice_line = None
        self.demplot = None
        self.plotmax = 0.0
        self.count = 0
        self.clicking = False
        self.last_click = 1
        self.slice_points = None
        self.dem_along_line = None
        self.max_line = None
        self.demimage = None
        self.slice_ticks = None
        self.logt = None

        self.red_temp = None
        self.grn_temp = None
        self.blu_temp = None
        self.colorbar = None
        self.the_algorithm = None
        self.crosshairs = []
        self.crosshairs_bezier = []
        self.demlines = []
        self.slice_points_interpolated = []
        self.slice_ticks_list = []
        self.dem_vertlines = []
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




    def displays(self, debug=False):
        # ui0 = widgets.HBox([self.rtemp,self.gtemp,self.btemp,self.sigma,]) if debug else widgets.HBox([])
        # ui05= widgets.HBox([self.slice_type, self.rng, self.tick_spacing])
        # ui11 = widgets.HBox([self.normalization, self.mouseover])
        # ui1 = widgets.HBox([self.algorithm])
        # ui15 = widgets.HBox([self.btn_draw_curve, self.btn_reset_lines])
        # ui = widgets.VBox([ui0,ui05, ui11, ui1, ui15])
        ui0 = widgets.HBox([self.rtemp,self.gtemp,self.btemp,self.sigma]) if debug else widgets.HBox([])
        #ui05= widgets.HBox([self.slice_type, self.rng, self.tick_spacing])
        #ui11 = widgets.HBox([self.normalization, self.mouseover])
        #ui1 = widgets.HBox([self.algorithm])
        ui15 = widgets.HBox([self.btn_draw_curve, self.btn_reset_lines,self.rng, self.algorithm, self.mouseover])
        ui = widgets.VBox([ui0, ui15])
        out = widgets.interactive_output(self.widgwrap, {'rtemp': self.rtemp, 'gtemp': self.gtemp, 'btemp': self.btemp, 'sigma': self.sigma,
                                        'algorithm': self.algorithm, 'rng': self.rng, 'slice_type': self.slice_type,
                                        "mouseover": self.mouseover, "spacing": self.tick_spacing, 'normalization': self.normalization})
        return ui, out

    def display(self, debug=False):
        ui, out = self.displays(debug)
        display(ui,out)

    def widgwrap(self, rtemp, gtemp, btemp, sigma, algorithm, rng, slice_type, mouseover, spacing, normalization):
        if self.fig is None:
            self.create_figure()
            self.init_figure( rtemp, gtemp, btemp, sigma, algorithm, rng=rng, slice_type=slice_type)
        else:
            self.count += 1
            self.update_figure( rtemp, gtemp, btemp, sigma, algorithm, rng=rng, slice_type=slice_type,
                               mouseover=mouseover, spacing=spacing, normalization=normalization)


    def create_figure(self):
        # print("Creating Dashboard")
        self.fontsize_prev = plt.rcParams.get('font.size')
        plt.rcParams.update({'font.size':self.fontsize})

        # Display the custom CSS in the notebook
        HTML(self.custom_css)
        self.fig = plt.figure(constrained_layout=True)
        #self.fig.set_size_inches(14, 10)
        self.fig.set_size_inches(self.xsize, self.ysize)

        spec = self.fig.add_gridspec(ncols=3, nrows=5, width_ratios=[0.1, 0.6, 0.6], height_ratios=[1.5, 1,1,1,1.5])
        self.ax1 = self.fig.add_subplot(spec[:, 0])
        self.ax2 = self.fig.add_subplot(spec[:-1, 1], projection=self.first.wcs)
        self.ax3 = self.fig.add_subplot(spec[0:2, 2])
        self.ax4 = self.fig.add_subplot(spec[-1, 1])
        self.ax5 = self.fig.add_subplot(spec[2:5, 2])
        spec.tight_layout


    def update_legend(self):
        # Hide the legend based on the condition
        self.legend = self.ax3.legend(loc='upper right', fontsize=12,bbox_to_anchor=(1, 1))
        self.fig.canvas.draw_idle()

    def init_dem_line(self, ix, iy):
        NC = self.count
        self.crosshairs.append(self.ax2.plot([ix], [iy], marker='+', color=f"C{NC}", markersize=25)[0])
        self.count += 1
        tt, dd = self.get_dem_at(ix, iy)
        themax = np.argmax(dd)
        the_max_temp = tt[themax]
        self.dem_vertlines.append(self.ax3.axvline(the_max_temp, color=f"C{NC}"))
        thelabel = f'Click {self.count} at [{ix:03}, {iy:03}, Max = {the_max_temp:0.2f}]' if self.count < 6 else None
        thelabel = f'{the_max_temp:0.2f}'
        self.demlines.append(self.ax3.plot(tt, dd, color=f"C{NC}", label=thelabel)[0])
        self.update_legend()

    def get_dem_at(self, ix, iy):
        [ptlogt, ptdem] = self.emc.compute_dem(ix, iy, logt=self.logt, algorithm=self.the_algorithm)
        return 10*ptlogt, ptdem/1.0e28

    def init_mouseover_line(self):
        NC = self.count
        self.crosshair_mouseover, = self.ax2.plot([], [], color='purple', marker='+', markersize=25)
        self.demplot_mouseover, = self.ax3.plot([],[], color='purple', ls="--", zorder=10000,  label=f"Offscreen")
        self.demplot_mouseover_vert = self.ax3.axhline(62, color='purple', ls="--", zorder=10000)
        self.update_legend()

    def init_figure(self, rtemp, gtemp, btemp, sigma, algorithm, gfac=1.0/2.2, plt_emmax=3.0e27, rng=[58, 68], slice_type="bezier", mouseover=True):
        # print("Initializing Dashboard")
        self.the_algorithm = algorithm
        self.the_slice_type = slice_type
        self.mouseover = mouseover
        self.logt = np.linspace(5.5, 7.5, 200)
        self.last_update_time = time.time()

        [synthchanlogts, synthchantresps] = lognormal_synthetic_channels([rtemp, gtemp, btemp], sigma)
        [cbcoll, cblogts, cbints, cbsigmas] = dem_color_table(synthchanlogts[0])

        [ilo, ihi, jlo, jhi] = [None]*4  # Update as needed
        synthdata = self.emc.synthesize_data(synthchanlogts, synthchantresps, ilo=ilo, ihi=ihi, jlo=jlo, jhi=jhi, algorithm=self.the_algorithm)
        self.demimage = np.stack([dat.data for dat in synthdata.data]).T
        colorbar_synthdata = cbcoll.synthesize_data(synthchanlogts, synthchantresps)
        clbimage = np.stack([dat.data for dat in colorbar_synthdata.data]).T
        self.demimg = self.ax2.imshow(((np.clip(self.demimage, 0, plt_emmax)/plt_emmax)**gfac).transpose((1, 0, 2)), interpolation="None")

        # [ptlogt, ptdem] = self.emc.compute_dem(0,0, algorithm=self.the_algorithm)
        self.init_mouseover_line()
        self.ax3.set_ylim(0, 1.1)
        self.ax3.set_xlim(*rng)
        self.ax5.set_ylim(*rng)

        self.update_legend()

        try:
            plt.suptitle(synthdata[0].meta['algorithm'] + ' inversion at ' + self.emc.data()[0].meta['date-obs'])
        except KeyError:
            plt.suptitle(synthdata[0].meta['ALGORITHM'] + ' inversion at ' + self.emc.data()[0].meta['date-obs'])

        self.ax2.set(title='RGB Composite DEM image')
        self.ax3.set(title='Diff Emission Measure', xlabel='Temperature (dB Kelvin)', ylabel='DEM (Mm [$10^9$ cm$^{-3}$]$^2$/dBK)')
        self.ax5.set(title='Diff Emission Measure', xlabel='Along the Line', ylabel='Temperature (dB Kelvin)')
        self.ax3.minorticks_on()
        self.ax4.set(title='RGB Composite DEM channel responses', xlabel='Temperature (dB Kelvin)')


        self.red_temp, = self.ax4.plot(10*synthchanlogts[0], synthchantresps[0], 'r')
        self.grn_temp, = self.ax4.plot(10*synthchanlogts[1], synthchantresps[1], 'g')
        self.blu_temp, = self.ax4.plot(10*synthchanlogts[2], synthchantresps[2], 'b')

        self.colorbar = self.ax1.imshow(((clbimage/np.max(clbimage))**gfac), aspect='auto', extent=[cbints[0], cbints[-1], 10*cblogts[0], 10*cblogts[-1]])
        self.ax1.set(title='Color Reference', ylabel='Temperature (dB Kelvin)', xlabel='Channel EM')

        self.init_interactivity()


    def init_buttons(self):
        self.n_control = 0
        self.btn_draw_curve = widgets.Button(description="Calculate Curve")
        self.btn_draw_curve.button_style = 'primary'
        self.btn_draw_curve.disabled = True
        self.btn_reset_lines = widgets.Button(description="Reset Lines")
        self.btn_reset_lines.disabled = True

    def init_interactivity(self):
        [nx,ny,nz] = self.demimage.shape

        def on_click(event):
            if event.inaxes == self.ax2:
                ix, iy = int(event.xdata), int(event.ydata)
                i = min(max(ix, 0), nx-1)
                j = min(max(iy, 0), ny-1)
                if self.drawing:
                    self.update_slice_curve(i, j)  # Function to draw/update the Bezier curve
                self.init_dem_line(i, j)
            self.update_legend()
        self.fig.canvas.mpl_connect('button_press_event', on_click)


        def on_mouseover(event):
            if self.mouseover:
                if self.demplot_mouseover is None:
                    self.init_mouseover_line()
                if event.inaxes == self.ax2:
                    ix, iy = int(event.xdata), int(event.ydata)
                    xlen, ylen, zlen = self.demimage.shape

                    if ix >= 0 and ix < xlen and iy >= 0 and iy < ylen:  # Check if ix and iy are within the bounds
                        self.crosshair_mouseover.set_data([ix],[iy])
                        self.demplot_mouseover.set_data(*self.get_dem_at(ix, iy))

                        themax = np.argmax(self.demplot_mouseover.get_ydata())
                        the_max_temp = self.demplot_mouseover.get_xdata()[themax]
                        self.demplot_mouseover.set_label(f"Mouse at [{ix}, {iy}] {the_max_temp:0.2f}")
                        self.demplot_mouseover.set_label(f"{the_max_temp:0.2f}")
                        self.demplot_mouseover_vert.remove()
                        self.demplot_mouseover_vert = self.ax3.axvline(the_max_temp, color='purple', ls="--", zorder=10000)
                        self.demplot_mouseover_vert.set_visible(True)

                        self.update_legend()
                else:
                    self.crosshair_mouseover.set_data([np.nan],[np.nan])
                    self.demplot_mouseover.set_data([np.nan],[np.nan])
                    self.demplot_mouseover.set_label("Offscreen")
                    self.demplot_mouseover_vert.set_visible(False)

                    self.update_legend()
        self.fig.canvas.mpl_connect('motion_notify_event', on_mouseover)

        def on_calc_curve_clicked(b):
            # print("Drawing curve")
            b.description = "Computing..."
            b.disabled = True
            try:
                self.update_slice_map()
                b.description += "done!"
                b.button_style = 'success'
            except IndexError as e:
                b.description += "Failed!"
                b.button_style = 'warning'
                raise e
                # with self.output:
                #     clear_output()

        self.btn_draw_curve.on_click(on_calc_curve_clicked)


        def on_reset_lines_clicked(b):
            for line in self.ax3.lines:
                line.remove()

            for crosshair in self.crosshairs:
                if crosshair is not None and crosshair.axes is not None:
                    crosshair.remove()
            self.crosshairs = []

            if self.slice_line is not None:
                self.slice_line.set_data([np.nan], [np.nan])

            for image in self.ax5.get_images():
                image.remove()

            if self.max_line is not None:
                if not isinstance(self.max_line, list):
                    print(self.max_line)
                    self.max_line.remove()
                    self.max_line = None

            for line in self.ax5.lines:
                line.remove()

            self.count = 0
            self.control_points = []
            self.ax3.set_prop_cycle(rcParams['axes.prop_cycle'])

            self.remove_slice_ticks()

            self.init_mouseover_line()

            self.update_legend()

            # Reset the "Draw Curve" button properties
            self.btn_draw_curve.description = "Calculate Curve"
            self.btn_draw_curve.button_style = 'primary'  # Default color
            self.btn_draw_curve.disabled = True
            b.disabled = True

        self.btn_reset_lines.on_click(on_reset_lines_clicked)

    def update_slice_map(self):
        if self.slice_points is None:
            self.remove_slice_ticks()
            return
        dems = [self.get_dem_at(int(np.round(i)), int(np.round(j))) for i, j in zip(self.slice_points[0], self.slice_points[1])]

        temperatures = dems[0][0]


        # print(self.the_normalization)


        if self.the_normalization == "area":
            func = np.sum
        elif self.the_normalization == "max":
            func = np.max
        elif self.the_normalization == "none":
            func = lambda x: 1

        the_map = np.stack([dem[1]/func(dem[1]) for dem in dems]).T

        max_line = [temperatures[np.argmax(the_map[:, i])] for i in range(len(dems))]
        self.max_line, = plt.step(max_line, 'r', where='post')

        self.dem_along_line = self.ax5.imshow(the_map, aspect='auto', extent=[0, len(dems), np.min(temperatures), np.max(temperatures)])

    def enable_buttons(self, both=True):
        self.btn_reset_lines.disabled = False
        self.btn_draw_curve.disabled = not both

    def update_slice_curve(self, i, j):
        if self.the_slice_type == "bezier":
            self.update_bezier_curve(i, j)

        elif self.the_slice_type == "spline":
            self.update_spline_curve(i, j)


    def add_control_point(self, i, j):
        self.control_points.append((i, j))
        self.crosshairs.append(self.ax2.plot([i], [j], marker='o', markeredgecolor=f"k", markersize=10)[0])
        self.n_control = len(self.control_points)
        self.enable_buttons(False)



    def update_bezier_curve(self, i, j):
        self.add_control_point(i, j)

        if len(self.control_points) < 2:
            return
        self.enable_buttons()

        def compute_bez_slice_points(control_points, num_points=252):
            n = len(control_points) - 1
            t = np.linspace(0, 1, num_points)
            slice_points = np.zeros((num_points, 2))
            for i in range(n + 1):
                binomial_coeff = scipy.special.comb(n, i)
                term = binomial_coeff * (t**i)[:, np.newaxis] * ((1 - t)**(n - i))[:, np.newaxis] * np.array(control_points[i])
                slice_points += term
            return slice_points

        self.slice_points = compute_bez_slice_points(self.control_points, num_points=252).T

        self.update_curve()


    def update_spline_curve(self, i, j):
        self.add_control_point(i, j)

        if len(self.control_points) < 4:
            return
        self.enable_buttons()

        def compute_spline_slice_points(control_points, num_points=252):
            # Get x and y coordinates of control points
            x_coords, y_coords = zip(*control_points)

            # Fit a cubic spline to the control points
            tck, u = splprep([x_coords, y_coords], s=0)

            # Evaluate the spline at a set of points to create the curve
            u_new = np.linspace(0, 1, num_points)
            spline_points = splev(u_new, tck)

            return spline_points

        self.slice_points = compute_spline_slice_points(self.control_points)

        self.update_curve()


    def interpolate_slice_points(self):
        xp0 = self.slice_points[0]
        yp0 = self.slice_points[1]

        # Compute the arc length along the curve
        L = np.zeros(xp0.shape)
        for i in range(1, len(xp0)):
            distance = np.sqrt((xp0[i] - xp0[i-1])**2 + (yp0[i] - yp0[i-1])**2)
            L[i] = L[i-1] + distance

        # Normalize the arc length to [0, 1]
        L /= L[-1]

        # Create a linearly spaced parameter 'L2' from 0 to 1
        L2 = np.linspace(0, 1, 252)

        # Interpolate 'xp0' and 'yp0' along the normalized arc length 'L' to obtain evenly spaced points 'xp2' and 'yp2'
        xp2 = UnivariateSpline(L, xp0)(L2)
        yp2 = UnivariateSpline(L, yp0)(L2)

        self.slice_points_interpolated = np.asarray([xp2, yp2])
        # 'xp2' and 'yp2' now contain the resampled points that are evenly spaced along the curve

        # Remove the old scatter plot
        self.remove_slice_ticks()
        self.make_slice_ticks()


    def remove_slice_ticks(self):
        if self.slice_ticks is not None:
            self.slice_ticks.remove()
            self.slice_ticks = None
            self.slice_ticks_list = []

        for line in self.ax5.lines:
            line.remove()

    def make_slice_ticks(self):
        for ii, pt in enumerate(self.slice_points_interpolated.T):
            if ii % self.tick_spacing_value == 0:
                self.slice_ticks_list.append(pt)

        self.slice_ticks_array = np.asarray(self.slice_ticks_list).T

        for ii, pt in enumerate(self.slice_ticks_list):
            self.ax5.axvline(ii*self.tick_spacing.value, color="cyan", ls=":", zorder=1000)

    def update_curve(self):
        if self.slice_points is not None:
            self.interpolate_slice_points()
            # Clear the previous curve and draw a new one
            if self.slice_line is not None:
                self.slice_line.set_data(self.slice_points_interpolated[0], self.slice_points_interpolated[1])
            else:
                self.slice_line, = self.ax2.plot(self.slice_points_interpolated[0], self.slice_points_interpolated[1], 'r-')

            self.slice_ticks = self.ax2.scatter(self.slice_ticks_array[0], self.slice_ticks_array[1], marker='o', color='cyan', s=20, zorder=100)
            self.fig.canvas.draw_idle()


    def update_figure(self, rtemp, gtemp, btemp, sigma, algorithm, gfac=1.0/2.2, plt_emmax=3.0e27,
                      rng=[55, 75], slice_type=None, mouseover=True, spacing=50, normalization="max"):
        # Update the plots using the stored handles

        # print("Updating Figure")
        self.ax3.set_xlim(*rng)
        self.ax5.set_ylim(*rng)

        if algorithm is not None:
            self.the_algorithm = algorithm
        if slice_type is not None:
            self.the_slice_type = slice_type
        self.tick_spacing_value = spacing
        self.the_normalization = normalization
        self.mouseover = mouseover

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
            synthdata = self.emc.synthesize_data(synthchanlogts, synthchantresps, ilo=ilo, ihi=ihi, jlo=jlo, jhi=jhi, algorithm=self.the_algorithm)
            demimage = np.stack([dat.data for dat in synthdata.data]).T
            self.demimg.set_data(((np.clip(demimage, 0, plt_emmax)/plt_emmax)**gfac).transpose((1, 0, 2)))

        try:
            plt.suptitle(synthdata[0].meta['algorithm'] + ' inversion at ' + self.emc.data()[0].meta['date-obs'])
        except KeyError:
            plt.suptitle(synthdata[0].meta['ALGORITHM'] + ' inversion at ' + self.emc.data()[0].meta['date-obs'])

        self.update_curve()
        self.update_slice_map()

        self.fig.canvas.draw_idle()



