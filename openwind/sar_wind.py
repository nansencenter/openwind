#!/usr/bin/env python
# Name:		    openwind.py
# Purpose:      Calculate wind speed from SAR images and wind direction
# Authors:      Morten Wergeland Hansen, Knut-Frode Dagestad
# License:      This file is part of OPENWIND. You can redistribute it or
#               modify under the terms of GNU General Public License, v.3
#               http://www.gnu.org/licenses/gpl-3.0.html

import os
import argparse
import warnings

import numpy as np

try:
    import matplotlib.pyplot as plt
except:
    print 'WARNING: Matplotlib not available, cannot make plots'

from nansat import Nansat, Nansatmap
from model_wind import ModelWind
from cmod5n import cmod5n_inverse

import pdb

class SARWind(Nansat, object):
    '''
    A class for calculating wind speed from SAR images using CMOD
    '''

    def __init__(self, sar_image, winddir=None, pixelsize=500, eResampleAlg=1):
        '''
            Parameters
            -----------
            sar_image : string
                        The SAR image filename - should be the original file.
            winddir :   int/float (an arbitratry wind direction),
                        numpy array (an array of wind directions, same size as
                        the SAR data),
                        string (the name of a file with wind field information
                            which can be opened by nansat),
                        Nansat (a Nansat object with wind direction),
                        None

                        Auxiliary wind field information needed to calculate
                        SAR wind (must be or have wind direction)
        '''
        if isinstance(sar_image, str) or isinstance(sar_image, unicode):
            super(SARWind, self).__init__(sar_image)
        elif isinstance(sar_image, Nansat):
            raise TypeError('Use of Nansat objects as input to SARWind is' \
                    ' disabled in the master branch.')
            #warnings.warn('Using Nansat object to calculate wind. Note that' \
            #        ' any previous reprojection is repeated to' \
            #        ' maintain correct azimuth.')
            #super(SARWind, self).__init__(sar_image.fileName)
            #self.reproject(sar_image)

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

        self.set_auxiliary(winddir)
        self.calculate_wind(eResampleAlg=eResampleAlg)

    def set_look_direction(self):
        # OBS - this will only work on unprojected data.
        # Also, azimuth_up may give result switched by 180 degrees from the
        # platform or sensor azimuth, at least if given file is netcdf of
        # previously exported nansat data... Needs fix
        # This function should be replaced by a band added by the SAR mappers -
        # see https://github.com/nansencenter/nansat/issues/57
        warnings.warn('Function set_look_direction only works for unprojected' \
                ' data. It will work if the input SAR image points to the' \
                ' original source, but not if it is an exported netcdf.')
        if self.get_metadata()['ANTENNA_POINTING'] == 'RIGHT':
            look_direction = self.azimuth_up() + 90
        else:
            look_direction = self.azimuth_up() - 90
        self.add_band(array=look_direction, parameters={
                        'name': 'sar_look_direction',
                        'time': self.get_time(self.sigma0_bandNo),
                })

    def reproject(self, *args, **kwargs):
        # Placeholder for future reprojection function
        # overloading Nansat.reproject(), after calculating
        # SAR look direction and storing as NumPy array
        if not self.has_band('sar_look_direction'):
            self.set_look_direction()
        super(SARWind, self).reproject(*args, **kwargs)

    def get_auxiliary(self, eResampleAlg=1):
        '''
            Get auxiliary information (Nansat object with wind direction)
            needed to estimate wind speed from SAR
        '''

        # TODO:
        # - if several instances,
        #       choose the one closest in time to SAR image
        # - should check if wind object really covers SAR image
        # - check that surface (10 m) winds are chosen
        # - allow use of Nansat objects containing winddirection,
        #       and not the U and V components

        if not self.winddir:
            try:
                mw = ModelWind(self.SAR_image_time, domain=self,
                               eResampleAlg=eResampleAlg)
            except ValueError as e:
                warnings.warn(e.message)
                return None
        else:
            mw = ModelWind(wind=self.winddir, domain=self,
                           eResampleAlg=eResampleAlg)

        return mw

    def set_auxiliary(self, winddir):
        '''
            Change current auxiliary information for wind speed calculation

            Parameters
            -----------
            winddir :   int/float (an arbitratry wind direction),
                        numpy array (an array of wind directions, same size as
                        the SAR data),
                        string (the name of a file with wind field information
                            which can be opened by nansat),
                        Nansat (a Nansat object with wind direction),
                        None

                        Auxiliary wind field information needed to calculate
                        SAR wind (must be or have wind direction)
        '''
        # check input type
        self.winddir=winddir

    def calculate_wind(self, winddir=None, storeModelSpeed=False,
                       eResampleAlg=1):
        '''
            Calculate wind speed from SAR sigma0 in VV polarization
        '''

        if winddir!=None:
            self.set_auxiliary(winddir)

        if not isinstance(self.winddir, int) and not isinstance(self.winddir,
                float) and not isinstance(self.winddir, np.ndarray):
            aux = self.get_auxiliary(eResampleAlg=eResampleAlg)
            if not aux:
                print 'Did not calculate wind speeds'
                return None
            winddir_time = aux.time

            # Check time difference between SAR image and wind direction object
            timediff = self.SAR_image_time - winddir_time
            try:
                secondsDiff = timediff.total_seconds()
            except: # for < python2.7
                secondsDiff = (timediff.microseconds +
                               (timediff.seconds + timediff.days *
                                24 * 3600) * 10**6) / 10**6
            hoursDiff = np.abs(secondsDiff/3600.)
            print 'Time difference between SAR image and wind direction: ' \
                    + '%.2f' % hoursDiff + ' hours'
            print 'SAR image time: ' + str(self.SAR_image_time)
            print 'Wind dir time: ' + str(winddir_time)
            if hoursDiff > 3:
                print '#########################################'
                print 'WARNING: time difference exceeds 3 hours!'
                print '#########################################'

            wind_u_bandNo = aux._get_band_number({
                                'standard_name': 'eastward_wind',
                            })
            wind_v_bandNo = aux._get_band_number({
                                'standard_name': 'northward_wind',
                            })
            # Get wind direction
            u_array = aux[wind_u_bandNo]
            v_array = aux[wind_v_bandNo]
            # 0 degrees meaning wind from North, 90 degrees meaning wind from East
            winddirArray = np.degrees(
                    np.arctan2(-u_array, -v_array))
        else:
            # Constant wind direction is input
            print 'Using constant wind (from) direction: ' + str(self.winddir) + \
                    ' degrees clockwise from North'
            if not np.shape(self.winddir):
                winddirArray = np.ones(self.shape())*self.winddir
            else:
                winddirArray = self.winddir
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
            self.add_band(array=aux['windspeed'], parameters={
                            'wkv': 'wind_speed',
                            'name': 'model_windspeed',
                            'time': winddir_time,
            })

        # TODO: Replace U and V bands with pixelfunctions
        u = -windspeed*np.sin((180.0 - winddirArray)*np.pi/180.0)
        v = windspeed*np.cos((180.0 - winddirArray)*np.pi/180.0)
        self.add_band(array=u, parameters={
                            'wkv': 'eastward_wind',
                            'time': winddir_time,
        })
        self.add_band(array=v, parameters={
                            'wkv': 'northward_wind',
                            'time': winddir_time,
        })

    def plot(self, numVectorsX = 20, show=True, clim=[3,10], scale=None,
            windspeedBand='windspeed', winddirBand='winddirection',
            northUp_eastRight=True):
        ''' Basic plotting function showing CMOD wind speed
        overlaid vectors in SAR image projection

        parameters
        ----------
        clim : list
            Color limits of the image.
        scale : None or float
            Data units per arrow length unit, e.g., m/s per plot width;
            a smaller scale parameter makes the arrow longer.
            If None, a simple autoscaling algorithm is used,
            based on the average vector length and the number of vectors.
            The arrow length unit is given by the scale_units parameter.

        '''

        try:
            sar_windspeed = self[windspeedBand]
        except:
            raise ValueError('SAR wind has not been calculated, ' \
                'execute calculate_wind(winddir) before plotting.')

        winddirReductionFactor = np.round(
                self.vrt.dataset.RasterXSize/numVectorsX)

        Vgeo = np.cos(np.radians(self[winddirBand]+180)) # add 180 deg to get "to"-direction
        Ugeo = np.sin(np.radians(self[winddirBand]+180))

        # Transformation of the wind directions
        Usat = ( Vgeo*np.sin(np.radians(self.azimuth_up())) -
                    Ugeo*np.cos(np.radians(self.azimuth_up())) )
        Vsat = -( Vgeo*np.cos(np.radians(self.azimuth_up())) +
                    Ugeo*np.sin(np.radians(self.azimuth_up())) )

        # Flip images if wanted
        if northUp_eastRight:
            if self.get_metadata()['ORBIT_DIRECTION'].lower()=='descending':
                # NOTE: 
                #       - the origin of ASAR grids is at first measurement (in
                #       time) at near range
                #       - the origin of Radarsat-2 grids is at first
                #       measurement (in time) but at far range
                #       - the origin of Sentinel-1 grids is ?
                #       - the origin of Cosmo-Skymed grids is ?
                if os.path.basename(self.fileName)[:3]=='ASA':
                    sar_windspeed = np.fliplr(sar_windspeed)
                Usat = -np.fliplr(Usat)
                Vsat = np.fliplr(Vsat)
            else:
                sar_windspeed = np.flipud(sar_windspeed)
                Usat = np.flipud(Usat)
                Vsat = -np.flipud(Vsat)

        X, Y = np.meshgrid(range(0, self.vrt.dataset.RasterXSize,
                                    winddirReductionFactor),
                           range(0, self.vrt.dataset.RasterYSize,
                                    winddirReductionFactor))

        Ux = Usat[Y, X]*10
        Vx = Vsat[Y, X]*10
        plt.imshow(sar_windspeed)
        plt.clim(clim)
        cbar = plt.colorbar()
        plt.quiver(X, Y, Ux, Vx, angles='xy', scale = scale)
        plt.axis('off')
        if show:
            plt.show()
        return plt

    def save_wind_map_image(self, fileName, scale=None, numArrowsRange=10,
                            landmask=True,
                            colorbar=True, fontsize=6,
                            drawgrid=True, tight=True, title=None, **kwargs):

        nMap = Nansatmap(self, resolution='l',figsize=(5, 8))

        if 'vmin' in kwargs.keys():
            vmin = kwargs.pop('vmin')
        else:
            vmin = 2
        if 'vmax' in kwargs.keys():
            vmax = kwargs.pop('vmax')
        else:
            vmax = 20

        # use wind direction "to" for calculating u and v
        winddirection = np.mod(self['winddirection'] + 180, 360)
        windspeed = self['windspeed']

        # Plot windspeed as colored grid
        nMap.pcolormesh(windspeed, vmin=vmin, vmax=vmax)

        # apply landmask to windspeeds
        windspeed = np.ma.masked_where( self.watermask()[1]==2, windspeed )

        quiPixelSpacing = int(np.round(windspeed.shape[1]/numArrowsRange))

        # specify the number of quiver vectors
        numQuiX = int(self.vrt.dataset.RasterXSize / quiPixelSpacing)
        numQuiY = int(self.vrt.dataset.RasterYSize / quiPixelSpacing)

        if scale is None and np.max(windspeed) <= 10.0:
            scale = 200
        elif scale is None:
            scale = 300

        # Draw quivers
        Ux = np.sin(np.radians(winddirection)) * windspeed
        Vx = np.cos(np.radians(winddirection)) * windspeed
        nMap.quiver(Ux, Vx, scale=scale,
                    quivectors=(numQuiY, numQuiX))
        # NO HARDCODING - GRID SIZES CHANGE
        #,
        #            width=0.002, X=0.8, Y=0.1, U=10, label='10 m/sec',
        #            fontproperties={'size':8})

        if colorbar:
            nMap.add_colorbar(fontsize=fontsize)

        if drawgrid:
            nMap.drawgrid()

        if title is not None:
            plt.title(title, fontsize=10)

        if tight:
            nMap.fig.tight_layout()

        nMap.save(fileName, landmask=landmask, **kwargs)


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

