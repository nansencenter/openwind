#-------------------------------------------------------------------------------
# Name:		bayes_wind.py
# Purpose:      
#
# Author:       Morten Wergeland Hansen
# Modified:	Morten Wergeland Hansen
#
# Created:	21.11.2014
# Last modified:08.07.2015 11:25
# Copyright:    (c) NERSC
# License:      
#-------------------------------------------------------------------------------
import numpy as np
from scipy.interpolate import griddata
from scipy.ndimage.filters import uniform_filter
from scipy.signal import medfilt2d

from nansat.nansat import Nansat
from nansat.domain import Domain
from nansat.nsr import NSR

from openwind.cmod5n import cmod5n_forward
from openwind.cdop import cdop
from openwind.sar_wind import SARWind

from geospaas.utils import nansat_filename
from geospaas.catalog.models import Dataset

def cost_function(apriori, obs, err):
    ''' 
    Return the cost function of `obs` with uncertainty `err` defined in the
    domain of `apriori`
    '''
    # NOTE: This function should be written in c (probably using numpy
    # libraries) to improve execution speed
    if not np.isscalar(obs) or not np.isscalar(err):
        raise IOError('Given observation and its uncertainty must be scalars')
    return np.square(( apriori - obs )/err)

def window_stdev(arr, radius):
    ''' standard deviation on sliding windows

    Code extract from
    http://stackoverflow.com/questions/18419871/improving-code-efficiency-standard-deviation-on-sliding-windows
    '''
    c1 = uniform_filter(arr, radius*2, mode='constant', origin=-radius)
    c2 = uniform_filter(arr*arr, radius*2, mode='constant', origin=-radius)
    std = ((c2 - c1*c1)**.5)[:-radius*2+1,:-radius*2+1]

    # adjust result to get it on the same grid as input arr
    grid_x, grid_y = np.meshgrid( range(np.shape(arr)[1]),
            range(np.shape(arr)[0]) )
    egrid_x, egrid_y = np.meshgrid(
            np.linspace(0.5+radius/2, np.shape(arr)[1]+0.5-radius-1, np.shape(arr)[1]-radius-1),
            np.linspace(0.5+radius/2, np.shape(arr)[0]+0.5-radius-1, np.shape(arr)[0]-radius-1))
    epoints = np.array([np.ndarray.flatten(egrid_x),
        np.ndarray.flatten(egrid_y)]).transpose()
    evalues = np.ndarray.flatten(std)

    return griddata(epoints, evalues, (grid_x, grid_y), method='cubic')

def grid_based_uncertainty(arr, radius):
    incr = radius+1
    aa = np.zeros((arr.shape[0], arr.shape[1]+incr*2))
    aa[:,incr:-incr] = arr
    aa[:,0:incr] = arr[:,-incr:]
    aa[:,-incr:] = arr[:,0:incr]
    err = window_stdev(aa,radius)

    return err[:,incr:-incr]

def LL2XY(EPSG, lon, lat):
    point = ogr.Geometry(ogr.wkbPoint)
    point.AddPoint(lon, lat)
    inSpatialRef = osr.SpatialReference()
    inSpatialRef.ImportFromEPSG(4326)
    outSpatialRef = osr.SpatialReference()
    outSpatialRef.ImportFromEPSG(EPSG)
    coordTransform = osr.CoordinateTransformation(inSpatialRef,
                            outSpatialRef)
    point.Transform(coordTransform)
    return point.GetX(), point.GetY()

class BayesianWind(SARWind):
    WIND_SPEED_RANGE = np.linspace(-20,20,81)
    MODEL_ERR = 1.5 # alternatively using method grid_based_uncertainty
    DOPPLER_ERR = 5
    # s0 error: use variance within grid cell (going from full res to, e.g.,
    # 500m) - this should not be constant..
    S0_ERR_FAC = 0.078 # Portabella et al. (1998, 2002)
    RESAMPLE_ALG = 1

    def __init__(self, filename, doppler_dataset=None, *args, **kwargs):

        super(BayesianWind, self).__init__(filename, *args, **kwargs)

        [u_apriori, v_apriori] = np.meshgrid(self.WIND_SPEED_RANGE,
                self.WIND_SPEED_RANGE)
        # 0 degrees is wind from North
        direction_apriori = (180./np.pi*np.arctan2(u_apriori, v_apriori) + 180.) % 360
        speed_apriori = np.sqrt(np.square(u_apriori) + np.square(v_apriori))

        # Get Nansat object of the model wind field
        #model_wind = self.get_source_wind(reprojected=False) # where did this
        # function go? anyway, the below should be equally fine..
        model_wind = Nansat(self.get_metadata('WIND_DIRECTION_SOURCE'))

        if doppler_dataset:
            nc_uris = doppler_dataset.dataseturi_set.filter(uri__endswith='.nc')
            fdg_tot = np.ones((len(nc_uris),self.shape()[0],self.shape()[1])) * np.nan
            fdg_err = np.ones((len(nc_uris),self.shape()[0],self.shape()[1])) * np.nan
            for ii, nc_uri in enumerate(nc_uris):
                fn = nansat_filename(nc_uri.uri)
                nn = Nansat(fn)
                try:
                    fdg = nn['fdg']
                    val = nn['valid_doppler']
                except ValueError:
                    val = nn['valid_scene']
                    dcp=nn['dop_coef_predicted']
                    dc=nn['dop_coef_observed']
                    rb=nn['range_bias_scene']
                    az=nn['azibias']
                    dc[dc>10000]=np.nan
                    dcp[dcp>10000]=np.nan
                    rb[rb>10000]=np.nan
                    az[az>10000]=np.nan
                    fdg = -(dc - dcp - rb - az)

                # Set invalid to nan...
                fdg[val==0]==np.nan
                # Estimate Doppler uncertainty
                fdg[fdg>100] = np.nan
                fdg[fdg<-100] = np.nan
                mask = np.isnan(fdg)
                fdg[mask] = np.interp(np.flatnonzero(mask), np.flatnonzero(~mask),
                    fdg[~mask])
                err_fdg = grid_based_uncertainty(fdg,2)
                err_fdg[err_fdg<self.DOPPLER_ERR]=self.DOPPLER_ERR
                nn.add_band(array=err_fdg, parameters={'name':'err_fdg'})
                nn.reproject(self, eResampleAlg=self.RESAMPLE_ALG, tps=True)
                try:
                    fdg = nn['fdg']
                    val = nn['valid_doppler']
                except ValueError:
                    val = nn['valid_scene']
                    dcp=nn['dop_coef_predicted']
                    dc=nn['dop_coef_observed']
                    rb=nn['range_bias_scene']
                    az=nn['azibias']
                    dc[dc>10000]=np.nan
                    dcp[dcp>10000]=np.nan
                    rb[rb>10000]=np.nan
                    az[az>10000]=np.nan
                    # TODO: check sign conventions for old doppler - this may be wrong for
                    # descending pass...
                    fdg = -(dc - dcp - rb - az)

                # Set invalid to nan...
                fdg[val==0]==np.nan
                fdg_tot[ii] = fdg
                fdg_err[ii] = nn['err_fdg']
            
            fdg_tot = np.nanmedian(fdg_tot,axis=0)
            # Do 2D spatial filtering... May perhaps be removed...
            fdg_tot = medfilt2d(fdg_tot)
            fdg_err = np.nanmedian(fdg_err,axis=0)
            fdg_err[fdg_tot<-100] = np.nan
            fdg_err[fdg_tot>100] = np.nan
            fdg_tot[fdg_tot<-100] = np.nan
            fdg_tot[fdg_tot>100] = np.nan
            #      - check cdop vs fdg
            fcdop = cdop(
                    self['model_windspeed'], 
                    np.mod(np.abs(self['look_direction']-self['winddirection']), 360),
                    self['incidence_angle'],
                    self.get_metadata(band_id=self.get_band_number({'standard_name':
                        'surface_backwards_scattering_coefficient_of_radar_wave'}),
                        key='polarization')
                )

        # Estimate sigma0 uncertainty
        s0 = self['sigma0_VV']
        err_s0 = self.S0_ERR_FAC*s0
        #    #mask = np.isnan(fdg)
        #    #fdg[mask] = np.interp(np.flatnonzero(mask), np.flatnonzero(~mask),
        #    #        fdg[~mask])
        #err_s0 = grid_based_uncertainty(s0,2)
        #import ipdb
        #ipdb.set_trace()
        #err_s0[err_s0<self.S0_ERR_FAC*s0] = self.S0_ERR_FAC*s0[err_s0<self.S0_ERR_FAC*s0]

        # Estimate model wind uncertainty (leads to adjustment near fronts)
        #model_px_resolution = int(np.round( 2 * model_wind.get_pixelsize_meters()[0] /
        #        self.get_pixelsize_meters()[0] ))
        uu = model_wind['U']
        vv = model_wind['V']
        err_u = grid_based_uncertainty(uu,2)
        err_v = grid_based_uncertainty(vv,2)

        model_wind.add_band(array=err_u, parameters={'name':'err_u'})
        model_wind.add_band(array=err_v, parameters={'name':'err_v'})

        # Reproject to SAR image
        model_wind.reproject(self, eResampleAlg=self.RESAMPLE_ALG, tps=True)

        # Get uncertainties in model wind  
        err_u = model_wind['err_u']
        # Without the below, uncertainties are lower in uniform areas - this
        # should be quite reasonable...
        #err_u[err_u<self.MODEL_ERR] = self.MODEL_ERR
        err_v = model_wind['err_v']
        #err_v[err_v<self.MODEL_ERR] = self.MODEL_ERR

        # Assign shape of SAR image to variable imshape
        imshape = self.shape()

        # Initialize result arrays
        ub_modcmod = np.ones(imshape)
        vb_modcmod = np.ones(imshape)
        ub_all = np.ones(imshape)
        vb_all = np.ones(imshape)

        self.has_doppler = np.zeros(imshape)

        model_u = model_wind['U']
        model_v = model_wind['V']
        sar_look = self[self._get_band_number({'standard_name':
                'sensor_azimuth_angle'})]
        inci = self['incidence_angle']
        pol = self.get_metadata(
            band_id=self.get_band_number({
                'standard_name': 'surface_backwards_scattering_coefficient_of_radar_wave'
                }), 
            key='polarization')
        # TODO: Try to increase speed using numba or cython..
        print 'Applying Bayesian on one-by-one pixel'
        for i in range(imshape[0]):
            print 'Row %d of %d'%(i+1,imshape[0])
            for j in range(imshape[1]):
                # There seems to be a problem with Radarsat-2 incidence angles
                # after resize (nan-values and erroneous resampling)
                # THIS IS NOT YET IN ANY GITHUB ISSUES...
                if np.isnan(inci[i,j]) or inci[i,j]<0 or s0[i,j]==0.0:
                    ub_modcmod[i,j] = np.nan
                    vb_modcmod[i,j] = np.nan
                    ub_all[i,j] = np.nan
                    vb_all[i,j] = np.nan
                    continue

                # Calculate model cost functions
                cost_model_u = cost_function(u_apriori, model_u[i,j], err_u[i,j])
                cost_model_v = cost_function(v_apriori, model_v[i,j], err_v[i,j])

                # Calculate sigma0 cost function
                cmod_s0 = cmod5n_forward(speed_apriori,
                        direction_apriori-sar_look[i,j],
                        np.ones(np.shape(speed_apriori))*inci[i,j])
                cost_sigma0 = cost_function( cmod_s0, s0[i,j], err_s0[i,j] )

                cost = cost_model_v + cost_model_u + cost_sigma0
                ind_min = np.where(cost==np.min(cost,axis=None))
                ub_modcmod[i,j] = u_apriori[ind_min]
                vb_modcmod[i,j] = v_apriori[ind_min]

                speedij = np.sqrt( np.square(ub_modcmod[i,j]) 
                        + np.square(vb_modcmod[i,j]) )
                # The Doppler will not be good for low wind speeds when the 
                # NRCS is small - this is a simple way to filter out that data.
                # TODO: improve method to mask pixels in sardoppler.py
                if speedij < 3:
                    fdg_tot[i,j] = np.nan

                if (doppler_dataset and 
                        not np.isnan(fdg_tot[i,j]) and 
                        fdg_err[i,j]!=0 and
                        not np.isnan(fdg_err[i,j])):
                    # Calculate Doppler cost function
                    self.has_doppler[i,j] = 1
                    cdop_fdg = cdop(speed_apriori, 
                        np.mod(np.abs(sar_look[i,j]-direction_apriori), 360),
                        np.ones(np.shape(speed_apriori))*inci[i,j],
                        pol
                    )
                    cost_doppler = cost_function(cdop_fdg, fdg_tot[i,j], fdg_err[i,j])
                    cost += cost_doppler

                    ind_min = np.where(cost==np.min(cost,axis=None))
                    ub_all[i,j] = u_apriori[ind_min]
                    vb_all[i,j] = v_apriori[ind_min]
                else:
                    # Set result to that of model and cmod result
                    ub_all[i,j] = np.nan #ub_modcmod[i,j]
                    vb_all[i,j] = np.nan #vb_modcmod[i,j]

                # Should give uncertainties as well
                #self.rms_u[i,j] = err_u[i,j] + err_v[i,j] + ...

        self.add_band(
            array = np.sqrt(np.square(ub_modcmod) + np.square(vb_modcmod)),
            parameters={
                'wkv': 'wind_speed',
                'name':'bspeed_modcmod',
                'long_name': 'Bayesian wind speed using model and ' \
                        'cmod data'}
            )
        self.add_band(
            array = np.mod(180. + 180./np.pi*np.arctan2(ub_modcmod,
                vb_modcmod), 360),
            parameters = {
                'wkv': 'wind_from_direction',
                'name': 'bdir_modcmod',
                'long_name': 'Bayesian wind direction using model and ' \
                        'cmod data'}
            )
        if doppler_dataset:
            self.add_band(
                array = fdg_tot,
                parameters={
                    'wkv': 'surface_backwards_doppler_frequency_shift_of_radar_wave_due_to_surface_velocity',
                    }
                )
            self.add_band(
                array = fdg_err,
                parameters={
                    'name': 'fdg_err',
                    }
                )
            self.add_band(
                array = fcdop,
                parameters={
                    'wkv': 'surface_backwards_doppler_frequency_shift_of_radar_wave_due_to_wind_waves',
                    'note': 'based on model wind'
                    }
                )
            self.add_band(
                array = np.sqrt(np.square(ub_all) + np.square(vb_all)),
                parameters={
                    'wkv': 'wind_speed',
                    'name':'bspeed_all',
                    'long_name': 'Bayesian wind speed using all data'}
                )
            self.add_band(
                array = np.mod(180. + 180./np.pi*np.arctan2(ub_all,
                    vb_all), 360),
                parameters = {
                    'wkv': 'wind_from_direction',
                    'name': 'bdir_all',
                    'long_name': 'Bayesian wind direction using all data'}
                )

    def export(self, *args, **kwargs):
        bands = kwargs.pop('bands', None)
        if not bands:
            bands = [
                    self.get_band_number('sigma0_VV'), 
                    self.get_band_number('model_windspeed'),
                    self.get_band_number('winddirection'), 
                    self.get_band_number('windspeed'),
                    self.get_band_number('bspeed_modcmod'), 
                    self.get_band_number('bdir_modcmod'),
                ]
            if self.has_band('fdg'):
                bands.append(self.get_band_number('fdg'))
                bands.append(self.get_band_number('fdg_err'))
                bands.append(self.get_band_number('fww'))
                bands.append(self.get_band_number('bspeed_all'))
                bands.append(self.get_band_number('bdir_all'))
        # TODO: add name of original file to metadata
        super(SARWind, self).export(bands=bands, *args, **kwargs)
