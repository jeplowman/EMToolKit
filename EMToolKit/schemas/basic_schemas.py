import os, sys, time, pickle, resource, copy, numpy as np
from EMToolKit.schemas import operators, element_functions, element_grid, coord_grid, element_source_responses
from EMToolKit.schemas.operators import multi_instrument_linear_operator, sparse_d2_partial_matrix
from EMToolKit.schemas.operators import reg_operator_postfac_wrapper, single_instrument_linear_operator_separable
from EMToolKit.schemas.basic_transforms import generic_transform, trivialframe, fits_transform
from EMToolKit.schemas.element_functions import (nd_voigt_psf, bin_function, get_2d_cov, get_3d_cov,
                               nd_gaussian_psf, nd_powgaussian_psf, spike_function,
                               flattop_guassian_psf, spice_spectrograph_psf)
from EMToolKit.schemas.element_grid import detector_grid, source_grid
from EMToolKit.schemas.coord_grid import coord_grid
from EMToolKit.schemas.element_source_responses import element_source_responses as esr
from sunpy.map import Map

# A detector schema needs to contain the information needed to map from a given source
# (or at least a source as defined by basic_source) onto its own detector numbers,
# as well as the code to compute the mapping.
class basic_detector(object):

	def __init__(self,meta):
		self.meta = meta
		self.wcs = None # Not yet implemented
		self.ndim = 2
		# We're not very consistent here about how many axes the source has. In some places we use
		# more general code. However, it elsewhere assumes we have a separable plane-of-sky+temperature
		# response style of detector, like AIA or XRT.
		self.transform = meta.get('TRANSFORM', fits_transform(self.meta))
		self.shape = meta['SHAPE']
		self.scale = np.array([np.sum((self.transform.index2coord(1,0)-self.transform.index2coord(0,0))**2)**0.5,
								np.sum((self.transform.index2coord(1,0)-self.transform.index2coord(0,0))**2)**0.5])
		#self.scale = np.array([meta.get('cdelt'+str(i)) for i in range(1,meta.get('naxis')+1)])
		#self.crpix = np.array([meta.get('crpix'+str(i)) for i in range(1,meta.get('naxis')+1)])
		#self.crval = np.array([meta.get('crval'+str(i)) for i in range(1,meta.get('naxis')+1)])
		self.origin = self.transform.index2coord(0.0,0.0)
		self.psfangle = meta.get('PSFANGLE',0.0)
		self.psfsigmas = np.array([meta.get('PSFSZ1',0.5),meta.get('PSFSZ2',0.5)])*self.scale
		self.psfcov = get_2d_cov(self.psfsigmas,self.psfangle)
		self.ipsfcov = np.linalg.inv(self.psfcov)
		# This assumes that the transform is affine:
		self.fwdtransform = np.matmul(self.transform.tfmat,np.diag(self.scale))
		self.frame = trivialframe(np.arange(1,self.ndim+1).astype(str))
		self.coords = coord_grid(self.shape,self.origin,self.fwdtransform,self.frame)
		self.grid = detector_grid(self.coords, [self.ipsfcov], nd_gaussian_psf)
		self.logt, self.tresp = meta['LOGT'], meta['TRESP']
		self.transform = meta.get('TRANSFORM',generic_transform)
		self.fwdops = []

	def fwdop(self,source):
		index = len(self.fwdops)
		for index in range(0,len(self.fwdops)):
			if(source.is_same(self.fwdops[index]['SOURCE'])): break
		if(index==len(self.fwdops)):
			sresp = esr(source.grid,self.grid,self.transform)
			sresp *= np.prod(source.scale)/np.prod(self.scale)/np.median(np.sum(sresp,axis=0).A1)
			operator = single_instrument_linear_operator_separable(sresp, self.tresp, temps=self.logt, exptime=self.meta['EXPTIME'])
			self.fwdops.append({'OPERATOR':operator,'SOURCE':source})
		return self.fwdops[index]['OPERATOR']

class basic_source(object):
	# This is a boneheaded source model that's based on an
	# ndcube data sequence (really anything that behaves like
	# a sunpy map should do. It assumes the map
	# has meta with crpix, cdelt, crval, and crota keywords
	# and does a basic transform assuming a fixed observer
	# coordinate. This should be updated to use the astropy WCS!!
	# Ideally the forward operation should be part of what the
	# data object does, but perhaps we don't want to load that
	# on the data objects since they might not be expected
	# to work on the multi-instrument DEM paradigm?
	# This is also a plane-of-sky+DEM model, and the modelheader
	# cdelt, etc parameters are written as if it's a simple image
	def __init__(self,sequence,super_fac=1, logt=None):
		nc = len(sequence)

		def minmax(arg):
			return([np.min(arg),np.max(arg)])

		date = None
		# Need to find out the smallest pixel size in the sequence and the maximum extent of the data in it
		for i in range(0,nc):
			logt = sequence[i].meta['LOGT']
			tmini,tmaxi = minmax(logt)
			dti = np.min(np.abs(logt[1:]-logt[0:-1]))
			nxi, nyi = sequence[i].data.shape
			bft = fits_transform(sequence[i].meta)
			corner00, corner01 = bft.index2coord(0,0), bft.index2coord(0,nyi)
			corner10, corner11 = bft.index2coord(nxi,0), bft.index2coord(nxi,nyi)
			xmini,xmaxi = minmax([corner00[0],corner01[0],corner10[0],corner11[0]])
			ymini,ymaxi = minmax([corner00[1],corner01[1],corner10[1],corner11[1]])
			dxi, dyi = np.sum((bft.index2coord(1,0)-corner00)**2)**0.5, np.sum((bft.index2coord(0,1)-corner00)**2)**0.5
			if(i>0):
				dx, dy, dt = [np.min([dxi,dx]),np.min([dyi,dy]),np.min([dti,dt])]
				tmin,xmin,ymin = [np.min([tmin,tmini]),np.min([xmin,xmini]),np.min([ymin,ymini])]
				tmax,xmax,ymax = [np.max([tmax,tmaxi]),np.max([xmax,xmaxi]),np.max([ymax,ymaxi])]
			else: [dx,dy,tmin,tmax,xmin,xmax,ymin,ymax,dt] = [dxi,dyi,tmini,tmaxi,xmini,xmaxi,ymini,ymaxi,dti]
			if(date is None): date=sequence[i].meta.get('DATE-OBS',sequence[i].meta.get('DATE-AVG'))
		dx, dy = dx*super_fac, dy*super_fac
		nx,ny = np.ceil((xmax-xmin)/dx).astype(np.uint32)+1, np.ceil((ymax-ymin)/dy).astype(np.uint32)+1

		check = True
		if(logt is None):
			if(nc == 1): logt = sequence[0].meta['LOGT']
			else: # Figure out if all of the logts are the same
				for i in range(1,nc):
					check *= len(sequence[i].meta['LOGT']) == len(sequence[i-1].meta['LOGT'])
					if(check): check *= np.sum((sequence[i].meta['LOGT']-sequence[i-1].meta['LOGT'])**2) == 0
		if(check): logt, nt = sequence[0].meta['LOGT'], len(sequence[0].meta['LOGT'])
		else:
			nt = np.round((tmax-tmin)/dt).astype(np.uint32)+1
			logt = tmin+dt*np.arange(nt)

		self.logt = logt
		basislogt = np.linspace(np.min(logt),np.max(logt),2*(nt-1)+1)
		[self.logts,self.bases] = [np.tile(basislogt,[nt,1]),np.zeros([nt,len(basislogt)])]
		for i in range(0,nt):
			self.bases[i] = (basislogt == logt[i])
			if(i > 0): self.bases[i] += (basislogt-logt[i-1])*(basislogt < logt[i])*(basislogt > logt[i-1])/(logt[i]-logt[i-1])
			if(i < nt-1): self.bases[i] += (logt[i+1]-basislogt)*(basislogt < logt[i+1])*(basislogt > logt[i])/(logt[i+1]-logt[i])

		x0, y0 = xmin - 0.5*((nx-1)*dx-(xmax-xmin)), ymin - 0.5*((ny-1)*dy-(ymax-ymin))
		self.meta = {'CDELT1':dx, 'CDELT2':dy, 'CROTA':0.0, 'CRPIX1':1, 'CRPIX2':1, 'RSUN_REF':sequence[0].meta['RSUN_REF'],
					 'CRVAL1':x0, 'CRVAL2':y0, 'NAXIS1':nx, 'NAXIS2':ny, 'DSUN_OBS':sequence[0].meta['DSUN_OBS'],
					 'CUNIT1':sequence[0].meta['CUNIT1'], 'CUNIT2':sequence[0].meta['CUNIT2'], 'HGLN_OBS':sequence[0].meta['HGLN_OBS'],
					 'ctype1':sequence[0].meta['ctype1'], 'CTYPE2':sequence[0].meta['CTYPE2'], 'HGLT_OBS':sequence[0].meta['HGLT_OBS'],
					 'LOGT':logt, 'PARENT_WCS':sequence[0].wcs, 'DATE-OBS':date, 'DATE-AVG':date}
		self.shape, self.axes = [nt,nx,ny], [logt,x0+dx*np.arange(nx),y0+dy*np.arange(ny)]
		self.scale = np.array([dx,dy])
		self.spatial_shape = np.array([nx,ny])
		dummymap = Map(np.zeros(self.spatial_shape,dtype=bool),self.meta)
		self.wcs = dummymap.wcs # This is wrong and a placeholder, since the model may have different dimensions and scale. Need to build a WCS instead.
		self.transform = fits_transform(self.meta)
		self.origin = self.transform.index2coord(0.0,0.0)
		self.ndim_spatial = len(self.origin)
		self.fwdtransform = np.matmul(self.transform.tfmat,np.diag(self.scale))
		self.frame = trivialframe(np.arange(1,self.ndim_spatial+1).astype(str))
		self.coords = coord_grid(self.spatial_shape,self.origin,self.fwdtransform,self.frame)
		self.grid = source_grid(self.coords,None,bin_function)

	def is_same(self,src):
		check = 'meta' in dir(src)
		if(check):
			check = ((self.meta['CDELT1'] == src.meta.get('CDELT1'))*
					(self.meta['CDELT2'] == src.meta.get('CDELT2'))*
					(self.meta['CROTA'] == src.meta.get('CROTA'))*
					(self.meta['CRPIX1'] == src.meta.get('CRPIX1'))*
					(self.meta['CRPIX2'] == src.meta.get('CRPIX2'))*
					(self.meta['CRVAL1'] == src.meta.get('CRVAL1'))*
					(self.meta['CRVAL2'] == src.meta.get('CRVAL2'))*
					(self.meta['NAXIS1'] == src.meta.get('NAXIS1'))*
					(self.meta['NAXIS2'] == src.meta.get('NAXIS2'))*
					(len(self.meta['LOGT']) == len(src.meta.get('LOGT',[]))))
		return check
