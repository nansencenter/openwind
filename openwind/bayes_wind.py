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

from nansat.nansat import Nansat

from openwind.cmod5n import cmod5n_forward
from openwind.cdop import cdop
from openwind.sar_wind import SARWind

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

class BayesianWind(SARWind):
    wind_speed_range = np.linspace(-20,20,81)
    model_err = 1.5 # alternatively using method grid_based_uncertainty
    doppler_err = 5
    # s0 error: use variance within grid cell (going from full res to, e.g.,
    # 500m) - this should not be constant..
    s0_err_fac = 0.078 # Portabella et al. (1998, 2002)
    resample_alg = 1

    def __init__(self, filename, doppler_file='', *args, **kwargs):

        super(BayesianWind, self).__init__(filename, *args, **kwargs)

        [u_apriori, v_apriori] = np.meshgrid(self.wind_speed_range,
                self.wind_speed_range)
        direction_apriori = 180./np.pi*np.arctan2(u_apriori, v_apriori) # 0 is wind towards North
        speed_apriori = np.sqrt(np.square(u_apriori) + np.square(v_apriori))

        # Get Nansat object of the model wind field
        #model_wind = self.get_source_wind(reprojected=False) # where did this
        # function go? anyway, the below should be equally fine..
        model_wind = Nansat(self.get_metadata('WIND_DIRECTION_SOURCE'))

        if doppler_file:
            # Get Nansat object of the range Doppler shift
            dop = Nansat(doppler_file)
            # Estimate Doppler uncertainty
            fdg = dop['dop_coef_observed'] - dop['dop_coef_predicted'] - \
                    dop['range_bias_scene'] - dop['azibias']
            fdg[fdg>100] = np.nan
            fdg[fdg<-100] = np.nan
            mask = np.isnan(fdg)
            fdg[mask] = np.interp(np.flatnonzero(mask), np.flatnonzero(~mask),
                    fdg[~mask])
            err_fdg = grid_based_uncertainty(fdg,2)
            err_fdg[err_fdg<self.doppler_err]=self.doppler_err
            dop.add_band(array=err_fdg, parameters={'name':'err_fdg'})
            dop.reproject(self, eResampleAlg=self.resample_alg, tps=True)
            fdg = dop['dop_coef_observed'] - dop['dop_coef_predicted'] - \
                    dop['range_bias_scene'] - dop['azibias']
            err_fdg = dop['err_fdg']
            #fdg_err = dop['range_bias_std_scene'] - this is not the uncertainty...
        
        # Estimate sigma0 uncertainty
        s0 = self['sigma0_VV']
        err_s0 = self.s0_err_fac*s0
        #    #mask = np.isnan(fdg)
        #    #fdg[mask] = np.interp(np.flatnonzero(mask), np.flatnonzero(~mask),
        #    #        fdg[~mask])
        #err_s0 = grid_based_uncertainty(s0,2)
        #import ipdb
        #ipdb.set_trace()
        #err_s0[err_s0<self.s0_err_fac*s0] = self.s0_err_fac*s0[err_s0<self.s0_err_fac*s0]

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
        model_wind.reproject(self, eResampleAlg=self.resample_alg, tps=True)

        # Get uncertainties in model wind  
        err_u = model_wind['err_u']
        # Without the below, uncertainties are lower in uniform areas - this
        # should be quite reasonable...
        #err_u[err_u<self.model_err] = self.model_err
        err_v = model_wind['err_v']
        #err_v[err_v<self.model_err] = self.model_err

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

                if (doppler_file and 
                        fdg[i,j]>-100 and 
                        fdg[i,j]<100 and 
                        err_fdg[i,j]!=0 and
                        not np.isnan(err_fdg[i,j])):
                    # Calculate Doppler cost function
                    self.has_doppler[i,j] = 1
                    cdop_fdg = cdop(speed_apriori, 
                        sar_look[i,j]-direction_apriori,
                        np.ones(np.shape(speed_apriori))*inci[i,j], 'VV')
                    cost_doppler = cost_function(cdop_fdg, fdg[i,j], err_fdg[i,j])
                    cost += cost_doppler

                    ind_min = np.where(cost==np.min(cost,axis=None))
                    ub_all[i,j] = u_apriori[ind_min]
                    vb_all[i,j] = v_apriori[ind_min]


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
        if doppler_file:
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
