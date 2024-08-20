import copy, time, scipy.special
import numpy as np, matplotlib.pyplot as plt, ipywidgets as widgets
from IPython.display import display, HTML, clear_output
from scipy.interpolate import make_interp_spline, splprep, splev, UnivariateSpline
from matplotlib import rcParams

from astropy.coordinates import SkyCoord
from astropy.nddata import StdDevUncertainty
from sunpy.map import Map
from ndcube import NDCube, NDCubeSequence, NDCollection

import EMToolKit.EMToolKit as emtk
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

        rt0,gt0,bt0 = kwargs.get('rtemp',6.0), kwargs.get('gtemp',6.2), kwargs.get('btemp',6.4)
        sg0 = 0.15 #0.5*np.mean(np.sort([rt0,gt0,bt0])[1:]-np.sort([rt0,gt0,bt0])[0:-1])
        self.the_normalization = kwargs.get('normalization',"area")

        [nx,ny] = em_collection.collection[em_collection.collection['models'][0]][0].data.shape

        self.rtemp=widgets.FloatSlider(min=5, max=7, value=rt0, step=0.05, description='rtemp', continuous_update=False)
        self.gtemp=widgets.FloatSlider(min=5, max=7, value=gt0, step=0.05, description='gtemp', continuous_update=False)
        self.btemp=widgets.FloatSlider(min=5, max=7, value=bt0, step=0.05, description='btemp', continuous_update=False)
        self.sigma=widgets.FloatSlider(min=0.025, max=0.5, value=sg0, step=0.01, description='sigma', continuous_update=False)

        self.rng=widgets.FloatRangeSlider(min=55, max=75, value=(58, 68), step=0.5, description='PlotRange', continuous_update=False)
        self.algorithm=widgets.Dropdown(options=self.emc.collection['models'], description='algorithm', continuous_update=False)
        self.init_buttons()
        self.xsize,self.ysize=kwargs.get('xsize',15),kwargs.get('ysize',8)
        self.fontsize = kwargs.get('fontsize',10)
        self.uninitialized=True
        self.mouseover = widgets.Checkbox( value=True, description='mouseover')

        self.tick_spacing_value = 50
        self.control_points = []
        self.drawing = True

        self.slice_line = None
        self.count = 0
        self.slice_points = None
        self.slice_ticks = None

        self.crosshairs = []
        self.demlines = []
        self.slice_ticks_list = []

    def displays(self, debug=False):

        ui0 = widgets.HBox([self.rtemp,self.gtemp,self.btemp,self.sigma]) if debug else widgets.HBox([])
        ui15 = widgets.HBox([self.btn_reset_lines,self.rng, self.algorithm, self.mouseover])
        ui = widgets.VBox([ui0, ui15])

        out = widgets.interactive_output(self.widgwrap, {'rtemp': self.rtemp, 'gtemp': self.gtemp, 'btemp': self.btemp, 'sigma': self.sigma,
                                        'algorithm': self.algorithm, 'rng': self.rng,
                                        "mouseover": self.mouseover})

        return ui, out

    def display(self, debug=False):
        ui, out = self.displays(debug)
        display(ui,out)
        print("Click on the image to populate the dashboard")

    def widgwrap(self, rtemp, gtemp, btemp, sigma, algorithm, rng, mouseover):
        if self.uninitialized:
            self.create_figure()
            self.init_figure( rtemp, gtemp, btemp, sigma, algorithm, rng=rng)
            self.uninitialized=False
        else:
            self.count += 1
        self.update_figure( rtemp, gtemp, btemp, sigma, algorithm, rng=rng,
                            mouseover=mouseover)

    def create_figure(self):
        # print("Creating Dashboard")
        self.fontsize_prev = plt.rcParams.get('font.size')
        plt.rcParams.update({'font.size':self.fontsize})

        self.fig = plt.figure(figsize=[self.xsize,self.ysize])

        self.fig.canvas.toolbar_visible = False
        self.fig.canvas.header_visible = False
        self.fig.canvas.footer_visible = False

        spec = self.fig.add_gridspec(ncols=3, nrows=4, width_ratios=[1, 8, 6], height_ratios=[1,1,1,1], hspace=1)

        self.ax1 = self.fig.add_subplot(spec[:, 0]) #colorbar
        self.ax2 = self.fig.add_subplot(spec[:-1, 1], projection=self.first.wcs) # Color x-y Image
        self.ax3 = self.fig.add_subplot(spec[0:-2, 2]) # 1D DEM Curves
        self.ax4 = self.fig.add_subplot(spec[-1:, 1]) # Synthetic Channel Responses
        self.ax5 = self.fig.add_subplot(spec[-2:, 2]) # 2D DEM Image (length, temperature)
        spec.tight_layout(self.fig,pad=1.5,rect=(0.01,0,1,1))

    def init_dem_lineplot(self, ix, iy):
        NC = self.count
        self.crosshairs.append(self.ax2.plot([ix], [iy], marker='+', color=f"C{NC}", markersize=25)[0])
        self.count += 1
        tt, dd = self.get_dem_at(ix, iy)
        themax = np.argmax(dd)
        the_max_temp = tt[themax]
        thelabel = f'Click {self.count} at [{ix:03}, {iy:03}, Max = {the_max_temp:0.2f}]' if self.count < 6 else None
        thelabel = f'{the_max_temp:0.2f}'
        self.demlines.append(self.ax3.plot(tt, dd, color=f"C{NC}", label=thelabel)[0])

    def get_dem_at(self, ix, iy):

        [ptlogt, ptdem] = self.emc.compute_dem(ix, iy, logt=self.logt, algorithm=self.the_algorithm)
        # return np.ones_like(ptdem)*ix, np.ones_like(ptdem)*ix
        return 10*ptlogt, ptdem/1.0e28

    def init_mouseover_line(self):
        self.crosshair_mouseover, = self.ax2.plot([], [], color='purple', marker='+', markersize=25)
        self.demplot_mouseover, = self.ax3.plot([],[], color='purple', ls="--", zorder=10000,  label=f"Offscreen")

    def init_figure(self, rtemp, gtemp, btemp, sigma, algorithm, gfac=1.0/2.2, plt_emmax=3.0e27, rng=[58, 68], slice_type="spline", mouseover=True, width=None):

        self.the_algorithm = algorithm
        self.the_slice_type = slice_type
        self.mouseover = mouseover
        self.logt = np.linspace(5.5, 7.5, 200)

        [synthchanlogts, synthchantresps] = lognormal_synthetic_channels([rtemp, gtemp, btemp], sigma)
        [cbcoll, cblogts, cbints, cbsigmas] = dem_color_table(synthchanlogts[0])

        [ilo, ihi, jlo, jhi] = [None]*4  # Update as needed
        synthdata = self.emc.synthesize_data(synthchanlogts, synthchantresps, ilo=ilo, ihi=ihi, jlo=jlo, jhi=jhi, algorithm=self.the_algorithm)
        self.demimage = np.stack([dat.data for dat in synthdata.data]).T
        colorbar_synthdata = cbcoll.synthesize_data(synthchanlogts, synthchantresps)
        clbimage = np.stack([dat.data for dat in colorbar_synthdata.data]).T
        self.demimg = self.ax2.imshow(((np.clip(self.demimage, 0, plt_emmax)/plt_emmax)**gfac).transpose((1, 0, 2)), interpolation="None")

        self.init_mouseover_line()
        self.ax3.set_ylim(0, 1.1)
        self.ax3.set_xlim(*rng)
        self.ax5.set_ylim(*rng)

        try:
            alglabel = synthdata[0].meta['algorithm'] + ' at ' + self.emc.data()[0].meta['date-obs']
        except KeyError:
            alglabel = synthdata[0].meta['ALGORITHM'] + ' at ' + self.emc.data()[0].meta['date-obs']

        self.ax2.set(title='RGB Composite from '+alglabel)
        self.ax3.set(title='Diff Emission Measure', xlabel='Temperature (dB Kelvin)', ylabel='DEM (Mm [$10^9$ cm$^{-3}$]$^2$/dBK)')
        self.ax5.set(title='Diff Emission Measure', xlabel='Along the Line', ylabel='Temperature (dB Kelvin)')
        self.ax3.minorticks_on()
        self.ax4.set(title='RGB Composite DEM channel responses', xlabel='Temperature (dB Kelvin)')

        self.colorbar = self.ax1.imshow(((clbimage/np.max(clbimage))**gfac), aspect='auto', extent=[cbints[0], cbints[-1], 10*cblogts[0], 10*cblogts[-1]])
        self.ax1.set(title='Color\nReference', ylabel='Temperature (dB Kelvin)', xlabel='Channel EM')

        self.red_temp, = self.ax4.plot(10*synthchanlogts[0], 1*synthchantresps[0], 'r')
        self.grn_temp, = self.ax4.plot(10*synthchanlogts[1], 1*synthchantresps[1], 'g')
        self.blu_temp, = self.ax4.plot(10*synthchanlogts[2], 1*synthchantresps[2], 'b')

        self.init_interactivity()



    def init_interactivity(self):
        [nx,ny,nz] = self.demimage.shape

        def on_click(event):
            if event.inaxes == self.ax2:
                ix, iy = int(event.xdata), int(event.ydata)
                i = min(max(ix, 0), nx-1)
                j = min(max(iy, 0), ny-1)

                self.update_slice_curve(i,j)  # Function to draw/update the Bezier curve
                self.init_dem_lineplot(i,j)
                self.update_slice_map()

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
                        self.demplot_mouseover.set_data(*self.get_dem_at(ix, iy))  #THIS MIGHT NEED A TRANSPOSE
            elif(self.demplot_mouseover is not None):
                self.demplot_mouseover.remove()
                self.crosshair_mouseover.remove()
                self.demplot_mouseover = None
                self.crosshair_mouseover = None

        self.fig.canvas.mpl_connect('motion_notify_event', on_mouseover)

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

            for line in self.ax5.lines:
                line.remove()

            self.count = 0
            self.control_points = []
            self.ax3.set_prop_cycle(rcParams['axes.prop_cycle'])

            self.remove_slice_ticks()
            self.slice_points = None
            self.crosshair_mouseover.remove()
            self.init_mouseover_line()

            self.reset_buttons()
            b.disabled = True

        self.btn_reset_lines.on_click(on_reset_lines_clicked)

    def reset_buttons(self):
        # Reset the "Draw Curve" button properties
        self.btn_reset_lines.button_style = 'info'
        self.btn_reset_lines.disabled = True

    def init_buttons(self):
        self.n_control = 0
        self.btn_reset_lines = widgets.Button(description="Reset Lines")
        self.reset_buttons()

    def enable_buttons(self, both=True):
        self.btn_reset_lines.disabled = False

    def update_slice_map(self):
        if self.slice_points is None:
            self.remove_slice_ticks()
            return
        dems = [self.get_dem_at(int(np.round(i)), int(np.round(j))) for i, j in zip(self.slice_points[0], self.slice_points[1])]

        temperatures = dems[0][0]

        if self.the_normalization == "area":
            # print("Normalizing to the Sum")
            func = np.sum
        elif self.the_normalization == "max":
            func = np.max
            # print("Normalizing to the Max")
        elif self.the_normalization == "none":
            # print("No Normalization")
            func = lambda x: 1

        the_map = np.stack([dem[1]/func(dem[1]) for dem in dems]).T

        self.dem_along_line = self.ax5.imshow(the_map, aspect='auto', extent=[0, len(dems), np.min(temperatures), np.max(temperatures)])



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
            self.ax5.axvline(ii*self.tick_spacing_value, color="cyan", ls=":", zorder=1000)

    def update_curve(self):
        if self.slice_points is not None:
            self.interpolate_slice_points()
            # Clear the previous curve and draw a new one
            if self.slice_line is not None:
                self.slice_line.set_data(self.slice_points_interpolated[0], self.slice_points_interpolated[1])
            else:
                self.slice_line, = self.ax2.plot(self.slice_points_interpolated[0], self.slice_points_interpolated[1], 'r-')

            self.slice_ticks = self.ax2.scatter(self.slice_ticks_array[0], self.slice_ticks_array[1], marker='o', color='cyan', s=20, zorder=100)


    def update_figure(self, rtemp, gtemp, btemp, sigma, algorithm, gfac=1.0/2.2, plt_emmax=3.0e27,
                    rng=[55, 75], slice_type=None, mouseover=True, spacing=50, normalization="none", width=None):
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
            alglabel = synthdata[0].meta['algorithm'] + ' at ' + self.emc.data()[0].meta['date-obs']
        except KeyError:
            alglabel = synthdata[0].meta['ALGORITHM'] + ' at ' + self.emc.data()[0].meta['date-obs']

        self.ax2.set(title='RGB Composite from '+alglabel)
        self.update_curve()
        self.update_slice_map()


