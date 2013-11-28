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
    A class for retrieving wind speed from SAR images using CMOD
    '''

    def __init__(self, SAR_image, winddir='online', *args, **kwargs):

        # If SAR image is given as a Nansat object, 
        # we must regenerate the SARwind object from its fileName
        if isinstance(SAR_image, Nansat):
            inputSize = SAR_image.shape()
            SAR_image = SAR_image.fileName
        super(SARWind, self).__init__(SAR_image, *args, **kwargs)
        # If input SAR image was already resized, we must also 
        # resize the regenerated Nansat object
        if 'inputSize' in locals() and inputSize != self.shape():
            self.resize(width=inputSize[1], height=inputSize[0])
        # TODO: should issue error if input SAR image has been 
        # reprojected, as SAR look direction cannot be calculated


        #########################################################
        # Check that this is a SAR image
        #########################################################

        # Check that the SAR image has VV pol NRCS
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

        # Calculate wind directly, if requested by user
        if winddir is not None or \
            (kwargs.has_key('calculate_wind') and 
             kwargs.pop('calculate_wind') == True):
                self.calculate_wind(winddir)


    def reproject(self):
        # Placeholder for future reporjection function
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
            self.raw.dataset.SetMetadataItem('Model wind field',
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

            # Temporary workaround due to bug in Nansat (#43)
            winddir_tmp = Nansat(winddir.vrt.fileName)
            winddir_tmp.fileName = winddir.fileName

            # Bi-linear interpolation of wind onto SAR image domain
            winddir_tmp.reproject(self, eResampleAlg=1)

            # Get wind direction
            u_array = winddir_tmp[wind_u_bandNo]
            v_array = winddir_tmp[wind_v_bandNo]
            winddirArray = np.degrees(
                    np.arctan2(-u_array, -v_array)) # 0 from North, 90 from East
            windspeedModel = winddir_tmp[windspeed_bandNo]
     
        
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
        #
        windspeedCMOD[np.where(np.isnan(windspeedCMOD))] = np.nan
        windspeedCMOD[np.where(np.isinf(windspeedCMOD))] = np.nan

        # Make new wind object where wind is stored
        # Necessary, since it is presently not possible to add
        # bands to resized or reprojected Nansat objects,
        # though this also makes saving of netCDF files easier
        self.wind = Nansat(domain = self)

        # Add wind speed and direction as bands
        self.wind.add_band(array=windspeedCMOD, parameters={
                        'wkv': 'wind_speed',
                        'name': 'sar_windspeed',
                        'time': self.get_time(self.sigma0_bandNo),
                        'winddir_time': winddir_time
                })
        self.wind.add_band(array=winddirArray, parameters={
                            'wkv': 'wind_from_direction',
                            'name': 'model_winddir',
                            'time': winddir_time
                })
        if isinstance(winddir, Nansat):
            self.wind.add_band(array=windspeedModel, parameters={
                            'wkv': 'wind_speed',
                            'name': 'model_windspeed',
                            'time': self.get_time(self.sigma0_bandNo),
                            'winddir_time': winddir_time
                    })

        # Copy all metadata from SAR image to wind object
        self.wind.set_metadata(self.get_metadata())

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

        if not hasattr(self, 'wind'):
            raise ValueError('SAR wind has not been calculated, ' \
                'execute calculate_wind(winddir) before plotting.')
        winddirReductionFactor = np.round(
                self.vrt.dataset.RasterXSize/numVectorsX)
        # model_winddir is direction from which wind is blowing
        winddir_relative_up = 360 - self.wind['model_winddir'] + \
                                    self.azimuth_up()
        X, Y = np.meshgrid(range(0, self.vrt.dataset.RasterXSize, 
                                    winddirReductionFactor),
                           range(0, self.vrt.dataset.RasterYSize, 
                                    winddirReductionFactor))
        Ux = np.sin(np.radians(winddir_relative_up[Y, X]))
        Vx = np.cos(np.radians(winddir_relative_up[Y, X]))
        plt.imshow(self.wind['sar_windspeed'])
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
                ' (integer, 0 for wind from North, 90 for wind from East etc.)')
    parser.add_argument('-n', dest='netCDF', 
            help='Export numerical output to NetCDF file')
    parser.add_argument('-f', dest='figure_filename', 
            help='Save wind plot as figure (e.g. PNG or JPEG)')
    parser.add_argument('-r', dest='resize_factor', default=0.1, 
            help='Resize factor for SAR image (default = 0.1)', type=float)
    args = parser.parse_args()

    if args.figure_filename is None and args.netCDF is None:
        raise ValueError('Not much point in calculating SAR wind '\
            'if output is not saved as figure or netCDF...' \
            '\nTa deg en bolle!')

    # Read SAR image
    if args.SAR_filename[-3:] == '.N1':
        # Temporary workaround for ASAR due to Nansat bug #43
        SAR_image_tmp = Nansat(args.SAR_filename)
        SAR_image = Nansat(SAR_image_tmp.vrt.fileName)
    else:
        SAR_image = Nansat(args.SAR_filename)

    # Resize
    print 'Resizing SAR image with factor of ' + str(args.resize_factor)
    SAR_image.resize(factor=args.resize_factor)

    # Get wind direction
    try:
        winddir = int(args.winddir)
    except:
        winddir = args.winddir
    s = SARWind(SAR_image, winddir)

    # Save figure
    if args.figure_filename is not None:
        print 'Saving output as figure: ' + args.figure_filename
        plt = s.plot(show=False)
        plt.savefig(args.figure_filename, bbox_inches='tight', dpi=300)

    # Save as netCDF file
    if args.netCDF is not None:
        print 'Saving output to netCDF file: ' + args.netCDF
        s.wind.export(args.netCDF)
