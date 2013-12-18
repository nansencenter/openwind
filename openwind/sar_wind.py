#!/usr/bin/env python
# Name:		    openwind.py
# Purpose:      Calculate wind speed from SAR images and wind direction
# Authors:      Morten Wergeland Hansen, Knut-Frode Dagestad
# License:      This file is part of OPENWIND. You can redistribute it or
#               modify under the terms of GNU General Public License, v.3
#               http://www.gnu.org/licenses/gpl-3.0.html

import os
import argparse

import numpy as np

try:
    import matplotlib.pyplot as plt
except:
    print 'WARNING: Matplotlib not available, cannot make plots'

from nansat import Nansat
from cmod5n import cmod5n_inverse


class SARWind(Nansat, object):
    '''
    A class for calculating wind speed from SAR images using CMOD
    '''

    def __init__(self, SAR_filename, winddir='online', pixelsize=500):

        # Call constructor of superclass Nansat
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

        # Calculate wind directly, if requested by user
        if winddir is not None:
            self.calculate_wind(winddir)


    def reproject(self):
        # Placeholder for future reprojection function
        # overloading Nansat.reproject(), after calculating
        # SAR look direction and storing as NumPy array
        raise TypeError('Reprojection not yet implemented!')


    def calculate_wind(self, winddir='online'):

        #########################################################
        # Check wind direction input
        #########################################################
        # Try to fetch NCEP wind data online 
        #   if wind direction is not given
        if winddir == 'online':
            print 'Downloading online NCEP GFS wind...'
            import download_NCEP_GFS
            winddir = download_NCEP_GFS.download_ncep(
                            self.get_time()[self.sigma0_bandNo - 1])
            if winddir is None:
                raise ValueError('No NCEP file available for time of ' \
                        'SAR image. Can not calculate SAR wind speed.')
            downloaded_wind_file = winddir

        if winddir == 'archive':
            # Use a function to find wind files in local archive
            # To be implemented
            raise ValueError('Archive functionality not yet implemented.')

        # Filename is given as input
        if isinstance(winddir, str):
            winddir = Nansat(winddir)

        # A Nansat object is given as input, or read from filename
        if isinstance(winddir, Nansat):
            self.vrt.dataset.SetMetadataItem('Model wind field',
                winddir.mapper[7:])
            # Check if Nansat object contains wind direction
            try:
                windspeed_bandNo = winddir._get_band_number(
                                {'standard_name': 'wind_speed'})
                wind_u_bandNo = winddir._get_band_number(
                                {'standard_name': 'eastward_wind'})
                wind_v_bandNo = winddir._get_band_number(
                                {'standard_name': 'northward_wind'})
            except:
                raise TypeError(winddir.fileName + +
                    ' does not contain wind direction')

            # Check time difference between SAR image and wind direction object
            try:
                winddir_time = winddir.get_time(
                        wind_u_bandNo).replace(tzinfo=None)
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
            except:
                winddir_time = None
                print '#############################################'
                print 'WARNING: time of wind direction is not given!'
                print '#############################################'


            # TODO:
            # - if several instances, 
            #       choose the one closest in time to SAR image
            # - check that surface (10 m) winds are chosen
            # - should check if wind object really covers SAR image
            # - allow use of Nansat objects containing winddirection,
            #       and not the U and V components

            # Reproject wind field onto SAR image

            # Bi-linear interpolation of wind onto SAR image domain
            winddir.reproject(self, eResampleAlg=1)

            # Get wind direction
            u_array = winddir[wind_u_bandNo]
            v_array = winddir[wind_v_bandNo]
            winddirArray = np.degrees(
                    np.arctan2(-u_array, -v_array)) # 0 from North, 90 from East
            windspeedModel = winddir[windspeed_bandNo]
     
        # Constant wind direction is input
        if isinstance(winddir, int):
            print 'Using constant wind (from) direction: ' + str(winddir) + \
                    ' degrees clockwise from North'
            winddirArray = np.ones(self.shape())*winddir
            winddir_time = None

        #########################################################
        # Calculate SAR wind with CMOD
        #########################################################
        print 'Calculating SAR wind with CMOD...'
        # Find look direction
        # Note: this method is only correct for unprojected SAR images!
        if self.get_metadata()['ANTENNA_POINTING'] == 'RIGHT':
            look_direction = self.azimuth_up() + 90
        else:
            look_direction = self.azimuth_up() - 90

        winddir_relative_look = np.mod(winddirArray - look_direction, 360)
        windspeedCMOD = cmod5n_inverse(self[self.sigma0_bandNo], 
            winddir_relative_look, self['incidence_angle'])

        windspeedCMOD[np.where(np.isnan(windspeedCMOD))] = np.nan
        windspeedCMOD[np.where(np.isinf(windspeedCMOD))] = np.nan

        # Add wind speed and direction as bands
        self.add_band(array=windspeedCMOD, parameters={
                        'wkv': 'wind_speed',
                        'name': 'sar_windspeed',
                        'time': self.get_time(self.sigma0_bandNo),
                        'winddir_time': winddir_time
                })
        self.add_band(array=winddirArray, parameters={
                            'wkv': 'wind_from_direction',
                            'name': 'model_winddir',
                            'time': winddir_time
                })
        if isinstance(winddir, Nansat):
            self.add_band(array=windspeedModel, parameters={
                            'wkv': 'wind_speed',
                            'name': 'model_windspeed',
                            'time': self.get_time(self.sigma0_bandNo),
                            'winddir_time': winddir_time
                    })

        # Delete NCEP if it was downloaded 
        try:
            os.remove(downloaded_wind_file)
            print 'Deleted downloaded wind file: ' + downloaded_wind_file
        except:
            pass

        # TODO: 
        # - add other CMOD versions than CMOD5


    def plot(self, numVectorsX = 20, show=True):
        ''' Basic plotting function showing CMOD wind speed
        overlaid vectors in SAR image projection'''

        try:
            sar_windspeed = self['sar_windspeed']
        except:
            raise ValueError('SAR wind has not been calculated, ' \
                'execute calculate_wind(winddir) before plotting.')

        winddirReductionFactor = np.round(
                self.vrt.dataset.RasterXSize/numVectorsX)
        # model_winddir is direction from which wind is blowing
        winddir_relative_up = 360 - self['model_winddir'] + \
                                    self.azimuth_up()
        X, Y = np.meshgrid(range(0, self.vrt.dataset.RasterXSize, 
                                    winddirReductionFactor),
                           range(0, self.vrt.dataset.RasterYSize, 
                                    winddirReductionFactor))
        Ux = np.sin(np.radians(winddir_relative_up[Y, X]))
        Vx = np.cos(np.radians(winddir_relative_up[Y, X]))
        plt.imshow(sar_windspeed)
        plt.clim([3, 10])
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
                ' (integer, 0 for wind from North, 90 for wind from East etc.) '
                'Use string "online" for automatic download of NCEP GFS winds')
    parser.add_argument('-n', dest='netCDF', 
            help='Export numerical output to NetCDF file')
    parser.add_argument('-f', dest='figure_filename', 
            help='Save wind plot as figure (e.g. PNG or JPEG)')
    parser.add_argument('-p', dest='pixelsize', default=500, 
            help='Pixel size for SAR wind calculation (default = 500 m)', 
                type=float)
    args = parser.parse_args()

    if args.figure_filename is None and args.netCDF is None:
        raise ValueError('Not much point in calculating SAR wind '\
            'if output is not saved as figure or netCDF...' \
            '\nTa deg en bolle!')

    # Read SAR image
    sw = SARWind(args.SAR_filename, pixelsize=args.pixelsize)

    # Get wind direction
    try:
        winddir = int(args.winddir)
    except:
        winddir = args.winddir
    # Calculate wind
    sw.calculate_wind(winddir)

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

