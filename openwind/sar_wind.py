# Name:		open_wind.py
# Purpose:      
# Authors:      Morten Wergeland Hansen (mortenwh)
# License:      This file is part of OPENWIND. You can redistribute it or
#               modify under the terms of GNU General Public License, v.3
#               http://www.gnu.org/licenses/gpl-3.0.html
from nansat import Nansat
import numpy as np
import scipy
import os

from model_wind import ModelWind
from cmod5n import cmod5n_inverse

import pdb

class SARWind(Nansat, object):
    '''
    A class for retrieving wind speed from SAR images using cmod
    '''

    pixel_size = 500. # meter

    def __init__(self, *args, **kwargs):
        if kwargs.has_key('winddir'):
            winddir = kwargs.pop('winddir')
        else:
            winddir = np.nan
        super(SARWind, self).__init__(*args,**kwargs)

        # Check that this is a SAR acquisition that has a NanSat mapper and
        # that it has VV pol nrcs
        names = [self.bands()[k]['name'] for k in self.bands().keys()] 
        if not 'sigma0_VV' in names:
            raise TypeError(self.fileName + 
                    ' does not have SAR NRCS in VV polarization')
        
        if os.path.splitext(self.fileName)[1]=='.nc':
            return

        if 'PixelFunctionType' in self.get_metadata(bandID='sigma0_VV').keys():
            # copy sigma0_VV to numpy array and add this as a separate band to
            # avoid problem with resizing pixelfunction bands.
            # This is a temporary workaround until the nansat resize issue
            # nansencenter/nansat#41 is fixed
            s0 = self['sigma0_VV']
            self.add_band(array=s0, parameters={
                'wkv': 'surface_backwards_scattering_coefficient_of_radar_wave',
                'name': 's0_src4wind',
            })
        # workaround for incidence angle is needed in any case...
        inci = self['incidence_angle'] # we probably also have another problem with the FillValue
        self.add_band(array=inci, parameters={
            'wkv': 'angle_of_incidence',
            'name': 'inci_src4wind',
        })
        self.vrt = self.raw.copy()

        # Now reduce size of image 
        try:
            line_spacing = float(self.get_metadata()['LINE_SPACING'])
        except:
            # LINE_SPACING not defined in ASAR images
            line_spacing = 50
        self.resize(line_spacing/self.pixel_size) 
        # copy reduced vrt to raw to get correct (reduced) size of the
        # model_wind object created below
        self.raw = self.vrt.copy()

        if np.isnan(winddir):
            # Get model wind field
            model_wind = ModelWind(self)
            winddir = model_wind['winddirection']
            # Add metadata:
            self.raw.dataset.SetMetadataItem('Model wind field',
                    model_wind.mapper[7:])
        else:
            winddir = np.ones(np.shape(self[1]))*winddir
    
        # NOTE 1: The look direction is defined in the center of the
        # domain clockwise from north. For longer domains, especially at high
        # latitudes, the look direction may vary a lot over the domain, and
        # using the center angle will be a coarse approximation.
        look_direction = float(self.vrt.dataset.GetMetadataItem('SAR_center_look_direction'))
        winddir_relative = np.mod(winddir - look_direction, 360)
        windspeed = self.get_cmod_wind(winddir_relative)

        # Add wind speed and direction as bands
        self.add_band(array=windspeed, parameters={
                        'wkv': 'wind_speed',
                        'name': 'windspeed',
                    })
        self.add_band(array=winddir, parameters={
                        'wkv': 'wind_from_direction',
                        'name': 'winddirection',
                    })
        self.add_band(array=winddir_relative, parameters={
                        'wkv': 'wind_from_direction',
                        'name': 'winddir_look_relative',
                        'long_name': 'Wind direction relative to SAR look direction (0 is upwind)',
                    })

        # Winddir is defined as the direction from which the wind is blowing
        # positive clockwise from north. Switch angles and coordinate system
        # (v=x, u=-y) to get it mathematically correct:
        u = -windspeed*np.sin((180.0 - winddir)*np.pi/180.0)
        v = windspeed*np.cos((180.0 - winddir)*np.pi/180.0)

        self.add_band(array=u, parameters={
                        'wkv': 'eastward_wind',
                    })
        self.add_band(array=v, parameters={
                        'wkv': 'northward_wind',
                    })

        # copy raw VRT object to the current vrt
        self.vrt = self.raw.copy()

    def get_cmod_wind(self, winddir_relative):
        if 'PixelFunctionType' in self.get_metadata(bandID='sigma0_VV').keys():
            # Temporary workaround of nansencenter/nansat#41
            s0 = self['s0_src4wind']
        else:
            s0 = self['sigma0_VV']
        inci = self['inci_src4wind']
        # consider changing cmod to use only one line of incidence angles
        windspeedcmod5 = cmod5n_inverse( s0, winddir_relative,
                inci )
        windspeedcmod5[np.where(np.isnan(s0))] = np.nan
        windspeedcmod5[np.where(np.isinf(s0))] = np.nan

        return windspeedcmod5

