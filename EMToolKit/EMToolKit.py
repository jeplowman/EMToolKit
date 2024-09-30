import copy
import numpy as np
from sunpy.map import Map
from ndcube import NDCube, NDCubeSequence, NDCollection
from astropy.coordinates import SkyCoord
from astropy.nddata import StdDevUncertainty, UnknownUncertainty
from EMToolKit.schemas.basic_schemas import basic_detector, basic_source
import astropy
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