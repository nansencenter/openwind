#!/usr/bin/env python
# Name:		    openwind.py
# Purpose:      Calculate wind speed from SAR images and wind direction
# Authors:      Morten Wergeland Hansen, Knut-Frode Dagestad
# License:      This file is part of OPENWIND. You can redistribute it or
#               modify under the terms of GNU General Public License, v.3
#               http://www.gnu.org/licenses/gpl-3.0.html

import argparse
from datetime import datetime

import numpy as np

try:
    import matplotlib.pyplot as plt
    from matplotlib.cm import jet
except:
    print 'WARNING: Matplotlib not available, cannot make plots'

from nansat import Nansat, Domain
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


    def calculate_wind(self, winddir=None, storeModelSpeed=True):
        # Calculate wind speed from SAR sigma0 in VV polarization
        if winddir:
            self.winddir = winddir
        if not isinstance(self.winddir, int):
            if self.winddir is None or self.winddir == 'online':
                self.modelWind = ModelWind(time=self.SAR_image_time)
            else:
                if isinstance(self.winddir, ModelWind):
                    self.modelWind = self.winddir
                else:
                    self.modelWind = ModelWind(wind=self.winddir)
            winddir_time = self.modelWind.time 

            # Bi-linear interpolation onto SAR image
            self.modelWind.reproject(self, eResampleAlg=1)

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
            storeModelSpeed = False # No windsped available, if direction given as number

        # Calculate SAR wind with CMOD
        # TODO: 
        # - add other CMOD versions than CMOD5
        print 'Calculating SAR wind with CMOD...'
        startTime = datetime.now()

        windspeed = cmod5n_inverse(self[self.sigma0_bandNo], 
                        np.mod(winddirArray - self['SAR_look_direction'], 360), 
                        self['incidence_angle'])
        print 'Calculation time: ' + str(datetime.now() - startTime)

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


    def _get_masked_windspeed(self, landmask=True, icemask=True):
        try:
            sar_windspeed = self['windspeed']
        except:
            raise ValueError('SAR wind has not been calculated, ' \
                'execute calculate_wind(winddir) first.')

        sar_windspeed[sar_windspeed<0] = 0
        palette = jet

        if landmask:
            try: # Land mask
                sar_windspeed = np.ma.masked_where(
                                    self.watermask()[1]==2, sar_windspeed)
                palette.set_bad([.3, .3, .3], 1.0) # Land is masked (bad)
            except:
                print 'Land mask not available'
        
        if icemask:
            try: # Ice mask
                try: # first try local file 
                    ice = Nansat('metno_local_hires_seaice_' + 
                            self.SAR_image_time.strftime('%Y%m%d'), 
                            mapperName='metno_local_hires_seaice')
                except: # otherwise Thredds
                    ice = Nansat('metno_hires_seaice:' + 
                            self.SAR_image_time.strftime('%Y%m%d'))
                ice.reproject(self)
                iceBandNo = ice._get_band_number(
                    {'standard_name': 'sea_ice_area_fraction'})
                sar_windspeed[ice[iceBandNo]>0] = -1
                palette.set_under('w', 1.0) # Ice is 'under' (-1)
            except:
                print 'Ice mask not available'

        return sar_windspeed, palette

    def write_geotiff(self, filename, landmask=True, icemask=True):

        sar_windspeed, palette = self._get_masked_windspeed(landmask, icemask)

        nansat_geotiff = Nansat(array=sar_windspeed, domain=self,
                                parameters = {'name': 'masked_windspeed',
                                              'minmax': '0 20'})
                        
        nansat_geotiff.write_geotiffimage(filename)


    def plot(self, filename=None, numVectorsX = 18, show=True,
                landmask=True, icemask=True, flip=True):
        ''' Basic plotting function showing CMOD wind speed
        overlaid vectors in SAR image projection'''

        try:
            sar_windspeed, palette = self._get_masked_windspeed(landmask, icemask)
        except:
            raise ValueError('SAR wind has not been calculated, ' \
                'execute calculate_wind(winddir) before plotting.')

        winddirReductionFactor = np.round(
                self.vrt.dataset.RasterXSize/numVectorsX)
        # model_winddir is direction from which wind is blowing
        winddir_relative_up = 360 - self['winddirection'] + \
                                    self.azimuth_up()
        indX = range(0, self.vrt.dataset.RasterXSize, winddirReductionFactor)
        indY = range(0, self.vrt.dataset.RasterYSize, winddirReductionFactor)
        X, Y = np.meshgrid(indX, indY)
        try: # scaling of wind vector length, if model wind is available
            model_windspeed = self['model_windspeed']
            model_windspeed = model_windspeed[Y, X]
        except:
            model_windspeed = np.ones(X.shape)

        Ux = np.sin(np.radians(winddir_relative_up[Y, X]))*model_windspeed
        Vx = np.cos(np.radians(winddir_relative_up[Y, X]))*model_windspeed

        # Make sure North is up, and east is right
        if flip == True:
            lon, lat = self.get_corners()
            if lat[0] < lat[1]:
                sar_windspeed = np.flipud(sar_windspeed)
                Ux = -np.flipud(Ux)
                Vx = -np.flipud(Vx)

            if lon[0] > lon[2]:
                sar_windspeed = np.fliplr(sar_windspeed)
                Ux = np.fliplr(Ux)
                Vx = np.fliplr(Vx)

        # Plotting
        figSize = sar_windspeed.shape
        legendPixels = 60.0
        legendPadPixels = 5.0
        legendFraction = legendPixels/figSize[0]
        legendPadFraction = legendPadPixels/figSize[0]
        dpi=100.0

        fig = plt.figure()
        fig.set_size_inches((figSize[1]/dpi, (figSize[0]/dpi)*
                                (1+legendFraction+legendPadFraction)))
        ax = fig.add_axes([0,0,1,1+legendFraction])
        ax.set_axis_off()
        plt.imshow(sar_windspeed, cmap=palette, interpolation='nearest')
        plt.clim([0, 20])
        cbar = plt.colorbar(orientation='horizontal', shrink=.85,
                     fraction=legendFraction, pad=legendPadFraction)
        cbar.ax.set_ylabel('[m/s]', rotation=0)
        cbar.ax.yaxis.set_label_position('right')
        ax.quiver(X, Y, Ux, Vx, angles='xy', width=0.004,
                    scale=numVectorsX*10, scale_units='width',
                    color=[.0, .0, .0], headaxislength=4)
        if show:
            fig.show()
        if filename is not None:
            fig.savefig(filename, pad_inches=0, dpi=dpi)
        return fig


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
    sw = SARWind(args.SAR_filename, winddir, pixelsize=args.pixelsize)

    # Save figure
    if args.figure_filename is not None:
        print 'Saving output as figure: ' + args.figure_filename
        plt = sw.plot(filename=args.figure_filename, show=False)

    # Save as netCDF file
    if args.netCDF is not None:
        print 'NetCDF export temporarily disabled'
        print 'Waiting for Nansat #47:'
        print 'https://github.com/nansencenter/nansat/issues/47'
        #print 'Saving output to netCDF file: ' + args.netCDF
        #sw.export_wind(args.netCDF)
