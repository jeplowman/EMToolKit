import copy
import numpy as np
from sunpy.map import Map
from ndcube import NDCube, NDCubeSequence, NDCollection
from astropy.coordinates import SkyCoord
from astropy.nddata import StdDevUncertainty, UnknownUncertainty
from astropy.nddata import StdDevUncertainty, UnknownUncertainty
from EMToolKit.schemas.basic_schemas import basic_detector, basic_source
import astropy
from scipy.interpolate import interp1d
from scipy.integrate import trapezoid

from scipy.interpolate import interp1d
from scipy.integrate import trapezoid


def em_data(maps, errs, logts, tresps, channels=None):
    cubes = []
    for i in range(len(maps)):
        mapi = copy.copy(maps[i])
        if 'CHANNEL' not in mapi.meta:
            if channels is None:
                mapi.meta['CHANNEL'] = mapi.meta['DETECTOR'] + mapi.meta['WAVE_STR']
            else:
                mapi.meta['CHANNEL'] = channels[i]
        mapi.meta['LOGT'] = logts[i]
        mapi.meta['TRESP'] = tresps[i]
        mapi.meta['SHAPE'] = np.array(mapi.data.shape)
        if 'SCHEMA' not in mapi.meta:
            mapi.meta['SCHEMA'] = basic_detector(mapi.meta)
        cubes.append(NDCube(mapi, uncertainty=errs[i]))
    return NDCubeSequence(cubes)
def em_data(maps, errs, logts, tresps, channels=None):
    cubes = []
    for i in range(len(maps)):
        mapi = copy.copy(maps[i])
        if 'CHANNEL' not in mapi.meta:
            if channels is None:
                mapi.meta['CHANNEL'] = mapi.meta['DETECTOR'] + mapi.meta['WAVE_STR']
            else:
                mapi.meta['CHANNEL'] = channels[i]
        mapi.meta['LOGT'] = logts[i]
        mapi.meta['TRESP'] = tresps[i]
        mapi.meta['SHAPE'] = np.array(mapi.data.shape)
        if 'SCHEMA' not in mapi.meta:
            mapi.meta['SCHEMA'] = basic_detector(mapi.meta)
        cubes.append(NDCube(mapi, uncertainty=errs[i]))
    return NDCubeSequence(cubes)


def dem_model(coeffs, logts, bases, coord_info, algorithm, wrapper, meta=None, wrapargs={}):
    nd = len(coeffs)
    dem_sequence = []
    if isinstance(coord_info, astropy.wcs.wcs.WCS):
        wcs = coord_info
    else:
        meta = coord_info
        if not isinstance(coord_info, dict):
            print('Warning in EMToolKit dem_model: coord_info is not a wcs or dict')
        wcs = meta.get('wcs', None)

    if meta is None:
        if wcs is not None:
            meta = dict(wcs.to_header())
        else:
            print('Warning in EMToolKit dem_model: need wcs or image meta')
    if 'LOGT' not in meta:
        meta['LOGT'] = logts[0]
    if 'SCHEMA' not in meta:
        meta['SCHEMA'] = basic_source([Map(coeffs[0], meta)])
    for i in range(nd):
        logt0 = np.median(logts[i][np.where(bases[i] == np.max(bases[i]))])
        metai = {
            'BASIS': bases[i],
            'LOGT': logts[i],
            'LOGT0': logt0,
            'ALGORITHM': algorithm,
            'WRAPPER': wrapper,
            'WRAPARGS': wrapargs,
        }
        dem_sequence.append(NDCube(coeffs[i], wcs=wcs, meta={**meta, **metai}))
    return NDCubeSequence(dem_sequence, meta={'ALGORITHM': algorithm, 'WRAPPER': wrapper, 'WRAPARGS': wrapargs})
    if meta is None:
        if wcs is not None:
            meta = dict(wcs.to_header())
        else:
            print('Warning in EMToolKit dem_model: need wcs or image meta')
    print(meta)
    if 'LOGT' not in meta:
        meta['LOGT'] = logts[0]
    if 'SCHEMA' not in meta:
        meta['SCHEMA'] = basic_source([Map(coeffs[0], meta)])
    for i in range(nd):
        logt0 = np.median(logts[i][np.where(bases[i] == np.max(bases[i]))])
        metai = {
            'BASIS': bases[i],
            'LOGT': logts[i],
            'LOGT0': logt0,
            'ALGORITHM': algorithm,
            'WRAPPER': wrapper,
            'WRAPARGS': wrapargs,
        }
        dem_sequence.append(NDCube(coeffs[i], wcs=wcs, meta={**meta, **metai}))
    return NDCubeSequence(dem_sequence, meta={'ALGORITHM': algorithm, 'WRAPPER': wrapper, 'WRAPARGS': wrapargs})


class em_collection:
    def __init__(self, datasequence):
        self.collection = NDCollection([("data", datasequence), ("models", [])])
        self.precomputed_interpolations = None

    def data(self):
        return self.collection['data']
    def __init__(self, datasequence):
        self.collection = NDCollection([("data", datasequence), ("models", [])])
        self.precomputed_interpolations = None

    def data(self):
        return self.collection['data']

    def add_model(self, modelsequence):
        pairs = [(k, self.collection.get(k)) for k in self.collection.keys()]

        # Safely retrieve the algorithm name, with a fallback to 'algorithm'
        algorithm_name = modelsequence.meta.get('ALGORITHM', modelsequence.meta.get('algorithm'))
        pairs.append((algorithm_name, modelsequence))

        # Clear and recreate the collection with the updated pairs
        self.collection.clear()
        self.collection = NDCollection(pairs)

        # Append the algorithm name to the models list
        if 'models' not in self.collection:
            self.collection['models'] = []
        self.collection['models'].append(algorithm_name)

    def precompute_interpolations(self): #TODO call this function at the correct time!
        """Precompute interpolation functions for all pixels."""
        if self.precomputed_interpolations is None:
            self.precomputed_interpolations = {}
        for algorithm in self.collection['models']:
            if algorithm not in self.precomputed_interpolations.keys():
                # print(f"Precomputing {algorithm}", end="\n")
                model = self.collection[algorithm]
                interpolations = []
                for component in model:
                    complogt = component.meta.get('LOGT', component.meta.get('logt'))
                    compbasis = component.meta.get('BASIS', component.meta.get('basis'))
                    interp_func = interp1d(complogt, compbasis, fill_value=0.0, bounds_error=False)
                    interpolations.append(interp_func)
                self.precomputed_interpolations[algorithm] = interpolations

    def compute_dem(self, i, j, logt=None, algorithm=None):
        if algorithm is None:
            algorithm = self.collection['models'][0]
        model = self.collection[algorithm]

        if logt is None:
            logt = model[0].meta.get('LOGT', model[0].meta.get('logt'))

        dem = np.zeros(logt.size)
        interpolations = self.precomputed_interpolations[algorithm]
        for component, interp_func in zip(model, interpolations):
            if j <= component.data.shape[1] and i <= component.data.shape[0]:
                dem += component.data[i, j] * interp_func(logt)

        return logt, dem

    def compute_dem_all(self, logt=None, algorithm=None):
        if algorithm is None:
            algorithm = self.collection['models'][0]
        model = self.collection[algorithm]

        if logt is None:
            logt = model[0].meta.get('LOGT', model[0].meta.get('logt'))

        dem = np.zeros([model[0].data.shape[0], model[0].data.shape[1], logt.size])
        interpolations = self.precomputed_interpolations[algorithm]

        for component, interp_func in zip(model, interpolations):
            dem += np.expand_dims(component.data, -1) * interp_func(logt)

        return logt, dem

    def synthesize_data(self, logts, tresps, algorithm=None, channels=None, ilo=0, ihi=-1, jlo=0, jhi=-1, meta=None):
        if logts[0].size == 1:
            [logts, tresps] = [[logts], [tresps]]
        ndata = len(logts)
        if algorithm is None:
            algorithm = self.collection['models'][0]
        if channels is None:
            channels = ['SYNTHDATA' + str(i) for i in range(ndata)]
        model = self.collection[algorithm]
        [synthmaps, syntherrs] = [[], []]
        self.precompute_interpolations()

        for i in range(ndata):
            synthdata = np.zeros(model[0].data[ilo:ihi, jlo:jhi].shape)
            syntherrs.append(UnknownUncertainty(np.zeros(model[0].data.shape) - 1))
            synthmap = copy.deepcopy(model[0])[ilo:ihi, jlo:jhi]
            synthmap.meta['NAXIS'] = 2
            synthmap.meta['ALGORITHM'] = algorithm
            synthmap.meta['CHANNEL'] = channels[i]
            datainterp = interp1d(logts[i], tresps[i], fill_value=0.0, bounds_error=False)
            for ind, component in enumerate(model):
                basisinterp = self.precomputed_interpolations[algorithm][ind] #[model.index(component)]
                logt = np.unique(np.hstack([component.meta['LOGT'], logts[i]]))
                coupling = trapezoid(datainterp(logt) * basisinterp(logt), x=logt)
                synthdata += coupling * component.data[ilo:ihi, jlo:jhi]
            synthmap.data[:] = synthdata[:]
            synthmaps.append(synthmap)
    def synthesize_data(self, logts, tresps, algorithm=None, channels=None, ilo=0, ihi=-1, jlo=0, jhi=-1, meta=None):
        if logts[0].size == 1:
            [logts, tresps] = [[logts], [tresps]]
        ndata = len(logts)
        if algorithm is None:
            algorithm = self.collection['models'][0]
        if channels is None:
            channels = ['SYNTHDATA' + str(i) for i in range(ndata)]
        model = self.collection[algorithm]
        [synthmaps, syntherrs] = [[], []]
        self.precompute_interpolations()

        for i in range(ndata):
            synthdata = np.zeros(model[0].data[ilo:ihi, jlo:jhi].shape)
            syntherrs.append(UnknownUncertainty(np.zeros(model[0].data.shape) - 1))
            synthmap = copy.deepcopy(model[0])[ilo:ihi, jlo:jhi]
            synthmap.meta['NAXIS'] = 2
            synthmap.meta['ALGORITHM'] = algorithm
            synthmap.meta['CHANNEL'] = channels[i]
            datainterp = interp1d(logts[i], tresps[i], fill_value=0.0, bounds_error=False)
            for ind, component in enumerate(model):
                basisinterp = self.precomputed_interpolations[algorithm][ind] #[model.index(component)]
                logt = np.unique(np.hstack([component.meta['LOGT'], logts[i]]))
                coupling = trapezoid(datainterp(logt) * basisinterp(logt), x=logt)
                synthdata += coupling * component.data[ilo:ihi, jlo:jhi]
            synthmap.data[:] = synthdata[:]
            synthmaps.append(synthmap)

        return em_data(synthmaps, syntherrs, logts, tresps, channels=channels)
        return em_data(synthmaps, syntherrs, logts, tresps, channels=channels)

    def synthesize_map(self, map, logt=None, tresp=None, algorithm=None, channel=None):
        if 'CHANNEL' not in map.meta:
            if channel is None:
                map.meta['CHANNEL'] = map.meta['DETECTOR'] + map.meta['WAVE_STR']
            else:
                map.meta['CHANNEL'] = channel
        if 'SCHEMA' not in map.meta:
            map.meta['LOGT'] = logt
            map.meta['TRESP'] = tresp
            map.meta['SHAPE'] = np.array(map.data.shape)
            map.meta['SCHEMA'] = basic_detector(map.meta)
        if algorithm is None:
            algorithm = self.collection['models'][0]
        source = self.collection[algorithm][0].meta['SCHEMA']
        coeffs = np.array([data.data for data in self.collection[algorithm].data]).flatten()
        output_map = copy.deepcopy(map)
        output_map.meta['CHANNEL'] += '_SYNTHETIC'
        output_map.data[:] = (((map.meta['SCHEMA']).fwdop(source)) * coeffs).reshape(map.data.shape)
        return output_map

# Beta DEM uncertainty estimation code. First, a lengthy preamble:
# The DEM forward problem is given by $M_j = \sum_j R_{ij}c_j$ with data $D_i$
# and associated uncertainties $\sigma_i$, where the DEM is determined by $c_j$
# e.g., as $\sum_j B_j(T)c_j$ and $R_{ij} = \int R_i(T)B_j(T)dT$ ($R_i$ being the
# continuous response functions provided for the observations in question): solve
# for $c_j$ minimizing $\chi^2=\sum_i\frac{(D_i-M_i)^2}{\sigma_i^2}$, while
# satisfying constraints like $\vec{c} > 0$ and minizing some regularization
# operator like $\sum{ij}c_i\Lambda_{ij}c_j$.

# The uncertainty problem is similar: 'Given a solution for which $M_i=D_i$,
# how much can the $c_j$ change (resulting in $M'_i = \sum_j R_{ij}(c_j + \Delta
# c_j)$) before the resulting $\chi'^2$ is a nominal threshold value (usually
# this is the number of data points, $n_\mathrm{data}$)?

# This question has multiple answers depending on interpretation. Do we mean
# 'how much can an individual $c_j$ change?', or all of them? Do they change
# in the same way, randomly/incoherently, or do they change in such a way as
# to counterbalance each other? The uncertainty will depend drastically on the
# version of this question -- in the last version, it can easily be inifite or
# undefined. If they all change in the same way, the resulting uncertainty is
# fairly small. If only one is changing, then the uncertainty is larger (in
# proportion to the number of $c_j$). Randomly falls between by the square
# root of the number of parameters. We use the 'one at a time' estimate since
# it is the most conservative of the ones that are straightforward to compute,
# as a nod to the potential ill-posedness of the fully inqualified question.
# Conveniently, this uncertainty estimate doesn't require knowledge of the
# solution and is independent of it! It does, however, depend on the definition
# of a basis element; we will assume temperature bins of width 0.1 in
# $\log_{10}(T)$ (1 Dex), see below.

# If $\Delta c_j$ is the only source of deviation from agreement (i.e, none of
# the other c in the coefficient vector are changing) then $\chi'^2=n_{data}$
# implies that $\sum_i \frac{(R_{ij}\Delta c_j)^2}{\sigma_i^2} = n_\mathrm{data}$

# Therefore the uncertainty estimate is just

# $\Delta c_j = \frac{\sqrt(n_\mathrm{data})}{\sqrt{\sum_i R_{ij}^2/\sigma_i^2}}$

# In terms of the original temperature response functions,
# $R_{ij} = \int R_i(T) B_j(T) dT$. For purposes of this uncertainty estimate,
# we take the basis functions to be top hats of width $\Delta logT=0.1$ at
# temperature $T_j$. Therefore $R_{ij} = R_i(T_j)\Delta logT$ and the
# uncertainty in a DEM component with that characteristic width is

# $\sigma E_j = \frac{\sqrt(n_\mathrm{data}}{\sqrt{\sum_i (R_i(T_j) \Delta lotT)^2/\sigma_i^2}}$

# The input temperatures to this may have any spacing, irrespective of
# $\Delta logT$. The output therefore is an independent uncertainty
# estimate for an effective temperature bin of width $\Delta logT$ at that
# temperature.
    def estimate_uncertainty(self, logt, dlogt=0.1, project_map=None):
        data = self.data()

        if(project_map is None): project_map = data[0]

        output_cube = np.zeros([project_map.data.shape[0],project_map.data.shape[1],len(logt)])

        for i in range(0,len(data)):
            dat_tresp = data[i].meta['TRESP']
            dat_logt = data[i].meta['LOGT']
            dati = copy.deepcopy(data[i])
            dati.data[:] = dati.uncertainty.array
            tresp = np.interp(logt,dat_logt,dat_tresp,right=0,left=0)
            data_reproj = dati.reproject_to(project_map.wcs)
            for j in range(0,len(logt)):
                output_cube[:,:,j] += (tresp[j]*dlogt/data_reproj.data)**2

        return np.sqrt(len(data)/output_cube)

