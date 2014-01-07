#!/usr/bin/env python
# Name:		    openwind.py
# Purpose:      Calculate wind speed from SAR images and wind direction
# Authors:      Morten Wergeland Hansen, Knut-Frode Dagestad
# License:      This file is part of OPENWIND. You can redistribute it or
#               modify under the terms of GNU General Public License, v.3
#               http://www.gnu.org/licenses/gpl-3.0.html

import argparse

import numpy as np

try:
    import matplotlib.pyplot as plt
except:
    print 'WARNING: Matplotlib not available, cannot make plots'

from nansat import Nansat
from model_wind import ModelWind
from cmod5n import cmod5n_inverse


class SARWind(Nansat, object):
    '''
    A class for calculating wind speed from SAR images using CMOD
    '''

    def __init__(self, SAR_filename, winddir=None, pixelsize=500):
        '''
            Parameters
            -----------
            SAR_filename : string
                SAR image filename (original, raw file)
            winddir : int, string, Nansat, None
                Auxiliary wind field information needed to calculate
                SAR wind (must be or have wind direction)
        '''
        super(SARWind, self).__init__(SAR_filename)

        # Check that this is a SAR image with VV pol NRCS
        try:
            self.sigma0_bandNo = self._get_band_number(
                            {'standard_name': 
            'surface_backwards_scattering_coefficient_of_radar_wave', 
                            'polarization': 'VV'})
        except:
            raise TypeError(self.fileName + 
                ' does not have SAR NRCS in VV polarization')

        self.SAR_image_time = self.get_time(
                self.sigma0_bandNo).replace(tzinfo=None)

        if pixelsize != 'fullres':
            print 'Resizing SAR image to ' + str(pixelsize) + ' m pixel size'
            self.resize(pixelsize=pixelsize)

        self.winddir = winddir
        if winddir is not None:
            self.calculate_wind()

    def set_look_direction(self):
        # This function should be replaced by a band added by the SAR mappers -
        # see https://github.com/nansencenter/nansat/issues/57
        if self.get_metadata()['ANTENNA_POINTING'] == 'RIGHT':
            look_direction = self.azimuth_up() + 90
        else:
            look_direction = self.azimuth_up() - 90
        self.add_band(array=look_direction, parameters={
                        'name': 'sar_look_direction',
                        'time': self.get_time(self.sigma0_bandNo),
                })

    def reproject(self, *args, **kwargs):
        # Overloading Nansat.reproject()
        # Calculate SAR look direction before reprojection, and store as
        # NumPy band to avoid problems due to modified domain orientation
        if not self.has_band('sar_look_direction'):
            self.set_look_direction()
        super(SARWind, self).reproject(*args, **kwargs)

    def calculate_wind(self, winddir=None, storeModelSpeed=False):
        # Calculate wind speed from SAR sigma0 in VV polarization
        if winddir:
            self.winddir = winddir

        if not isinstance(self.winddir, int):
            if self.winddir is None:
                self.modelWind = ModelWind(self.SAR_image_time, domain=self)
            else:
                self.modelWind = ModelWind(wind=self.winddir, domain=self)
            winddir_time = self.modelWind.time 

            # Check time difference between SAR image and wind direction object
            timediff = self.SAR_image_time - winddir_time
            hoursDiff = np.abs(timediff.total_seconds()/3600.)
            print 'Time difference between SAR image and wind direction: ' \
                    + '%.2f' % hoursDiff + ' hours'
            print 'SAR image time: ' + str(self.SAR_image_time)
            print 'Wind dir time: ' + str(winddir_time)
            if hoursDiff > 3:
                print '#########################################'
                print 'WARNING: time difference exceeds 3 hours!'
                print '#########################################'

            wind_u_bandNo = self.modelWind._get_band_number({
                                'standard_name': 'eastward_wind',
                            })
            wind_v_bandNo = self.modelWind._get_band_number({
                                'standard_name': 'northward_wind',
                            })
            # Get wind direction
            u_array = self.modelWind[wind_u_bandNo]
            v_array = self.modelWind[wind_v_bandNo]
            winddirArray = np.degrees(
                    np.arctan2(-u_array, -v_array)) # 0 from North, 90 from East
        else:
            # Constant wind direction is input
            print 'Using constant wind (from) direction: ' + str(self.winddir) + \
                    ' degrees clockwise from North'
            winddirArray = np.ones(self.shape())*self.winddir
            winddir_time = None

        # Calculate SAR wind with CMOD
        # TODO: 
        # - add other CMOD versions than CMOD5
        print 'Calculating SAR wind with CMOD...'
        if not self.has_band('sar_look_direction'):
            self.set_look_direction()

        windspeed = cmod5n_inverse(self[self.sigma0_bandNo], 
                            np.mod(winddirArray - self['sar_look_direction'], 360), 
                            self['incidence_angle'])

        windspeed[np.where(np.isnan(windspeed))] = np.nan
        windspeed[np.where(np.isinf(windspeed))] = np.nan

        # Add wind speed and direction as bands
        # TODO: make it possible to update existing bands... See
        # https://github.com/nansencenter/nansat/issues/58
        self.add_band(array=windspeed, parameters={
                        'wkv': 'wind_speed',
                        'name': 'windspeed',
                        'time': self.get_time(self.sigma0_bandNo),
                        'winddir_time': winddir_time
                })
        self.add_band(array=winddirArray, parameters={
                            'wkv': 'wind_from_direction',
                            'name': 'winddirection',
                            'time': winddir_time
                })

        if storeModelSpeed:
            self.add_band(array=self.modelWind['windspeed'], parameters={
                            'wkv': 'wind_speed',
                            'name': 'model_windspeed',
                            'time': winddir_time,
            })

        # TODO: Replace U and V bands with pixelfunctions
        u = -windspeed*np.sin((180.0 - winddirArray)*np.pi/180.0)
        v = windspeed*np.cos((180.0 - winddirArray)*np.pi/180.0)
        self.add_band(array=u, parameters={
                            'wkv': 'eastward_wind',
        })
        self.add_band(array=v, parameters={
                            'wkv': 'northward_wind',
        })

    def plot(self, numVectorsX = 20, show=True):
        ''' Basic plotting function showing CMOD wind speed
        overlaid vectors in SAR image projection'''

        try:
            sar_windspeed = self['windspeed']
        except:
            raise ValueError('SAR wind has not been calculated, ' \
                'execute calculate_wind(winddir) before plotting.')

        winddirReductionFactor = np.round(
                self.vrt.dataset.RasterXSize/numVectorsX)
        # model_winddir is direction from which wind is blowing
        winddir_relative_up = 360 - self['winddirection'] + \
                                    self.azimuth_up()
        X, Y = np.meshgrid(range(0, self.vrt.dataset.RasterXSize, 
                                    winddirReductionFactor),
                           range(0, self.vrt.dataset.RasterYSize, 
                                    winddirReductionFactor))
        Ux = np.sin(np.radians(winddir_relative_up[Y, X]))
        Vx = np.cos(np.radians(winddir_relative_up[Y, X]))
        plt.imshow(sar_windspeed)
        plt.clim([0, 18])
        cbar = plt.colorbar()
        plt.quiver(X, Y, Ux, Vx, angles='xy')
        plt.axis('off')
        if show:
            plt.show()
        return plt


###################################
#    If run from command line
###################################
if __name__ == '__main__':

    parser = argparse.ArgumentParser()
    parser.add_argument('-s', dest='SAR_filename', 
            required=True, help='SAR image filename')
    parser.add_argument('-w', dest='winddir', 
            default='online', help='Wind direction filename or constant '
                ' (integer, 0 for wind from North, 90 for wind from East etc.). '
                'Omit this argument for automatic download of NCEP GFS winds.')
    parser.add_argument('-n', dest='netCDF', 
            help='Export numerical output to NetCDF file')
    parser.add_argument('-f', dest='figure_filename', 
            help='Save wind plot as figure (e.g. PNG or JPEG)')
    parser.add_argument('-p', dest='pixelsize', default=500, 
            help='Pixel size for SAR wind calculation (default = 500 m)', 
                type=float)
    args = parser.parse_args()

    if args.figure_filename is None and args.netCDF is None:
        raise ValueError('Please add filename of processed figure (-f) or' \
                ' netcdf (-n)')

    # Get wind direction
    try:
        winddir = int(args.winddir)
    except:
        winddir = args.winddir

    # Read SAR image
    sw = SARWind(args.SAR_filename, pixelsize=args.pixelsize)

    # Save figure
    if args.figure_filename is not None:
        print 'Saving output as figure: ' + args.figure_filename
        plt = sw.plot(show=False)
        plt.savefig(args.figure_filename, bbox_inches='tight', dpi=300)

    # Save as netCDF file
    if args.netCDF is not None:
        print 'NetCDF export temporarily disabled'
        print 'Waiting for Nansat #47:'
        print 'https://github.com/nansencenter/nansat/issues/47'
        #print 'Saving output to netCDF file: ' + args.netCDF
        #sw.export_wind(args.netCDF)
