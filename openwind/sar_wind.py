#!/usr/bin/env python
# Name:		openwind.py
# Purpose:      Calculate wind speed from SAR images and wind direction
# Authors:      Morten Wergeland Hansen, Knut-Frode Dagestad
# License:      This file is part of OPENWIND. You can redistribute it or
#               modify under the terms of GNU General Public License, v.3
#               http://www.gnu.org/licenses/gpl-3.0.html

from __future__ import absolute_import

import argparse
import warnings
from datetime import datetime

import numpy as np

try:
    import matplotlib.pyplot as plt
    from matplotlib.cm import jet
except:
    print 'WARNING: Matplotlib not available, cannot make plots'

from nansat.nansat import Nansat, Domain, _import_mappers
from nansat.nansatmap import Nansatmap
from openwind.cmod5n import cmod5n_inverse

class SARWind(Nansat, object):
    '''
    A class for calculating wind speed from SAR images using CMOD
    '''

    def __init__(self, sar_image, wind_direction='ncep_wind_online',
                    pixelsize=500, eResampleAlg=1, *args, **kwargs):
        '''
            Parameters
            -----------
            sar_image : string or Nansat object
                        The SAR image as a filename or Nansat object
            wind_direction : int, numpy array, string, Nansat
                        Auxiliary wind field information needed to calculate
                        SAR wind (must be or have wind direction in degrees):

                        - constant wind direction (integer),
                        - array of wind directions, same size as the SAR data,
                        - the name of a Nansat compatible file containing
                          wind direction information
                        - name of a mapper with functionality to find a wind
                          file (online or on local disk) matching the SAR
                          image time [DEFAULT: 'ncep_wind_online']
                        - a Nansat object with wind direction.
            pixel_size : float or int
                        Grid pixel size in metres
            eResampleAlg : int
                        Resampling algorithm used for reprojecting wind field
                        to SAR image
                            -1 : Average,
                             0 : NearestNeighbour
                             1 : Bilinear (default),
                             2 : Cubic,
                             3 : CubicSpline,
                             4 : Lancoz

        '''
        if isinstance(sar_image, str) or isinstance(sar_image, unicode):
            super(SARWind, self).__init__(sar_image, *args, **kwargs)
        elif isinstance(sar_image, Nansat):
            super(SARWind, self).__init__(domain=sar_image, *args, **kwargs)
            self.vrt = sar_image.vrt
            self.mapper = sar_image.mapper

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
        self._set_wind_direction_source(wind_direction)

        if pixelsize != 'fullres':
            print 'Resizing SAR image to ' + str(pixelsize) + ' m pixel size'
            self.resize(pixelsize=pixelsize)

        # If this is a netcdf file with already calculated windspeed (e.g.
        # created as a SARWind object in order to use the plotting functions),
        # do not calculate wind over again
        if not self.has_band('windspeed'):
            self.calculate_wind(wind_direction,eResampleAlg=eResampleAlg)

    def _set_wind_direction_source(self, wind_direction):
        '''
            Adds wind direction source information to metadata

            Parameters
            -----------
            wind_direction : int, numpy array, string, Nansat
                        Auxiliary wind field information needed to calculate
                        SAR wind (must be or have wind direction in degrees):

                        - constant wind direction (integer),
                        - array of wind directions, same size as the SAR data,
                        - the name of a Nansat compatible file containing
                          wind direction information
                        - name of a mapper with functionality to find a wind
                          file (online or on local disk) matching the SAR
                          image time [DEFAULT: 'ncep_wind_online']
                        - a Nansat object with wind direction.
        '''

        if isinstance(wind_direction, str):
            import nansat.nansat
            mnames = [key.replace('mapper_','') for key in nansat.nansat.nansatMappers]
            if wind_direction in mnames:
                self.set_metadata('WIND_DIRECTION_SOURCE', 
                    wind_direction + datetime.strftime(
                                    self.SAR_image_time, ':%Y%m%d%H%M'))
            else:
                self.set_metadata('WIND_DIRECTION_SOURCE',
                    wind_direction)
        elif isinstance(wind_direction,Nansat):
            self.set_metadata('WIND_DIRECTION_SOURCE',
                    'Nansat instance of '+wind_direction.fileName)
        elif isinstance(wind_direction, int):
            self.set_metadata('WIND_DIRECTION_SOURCE',
                    'Set by user ('+str(wind_direction)+')')
        elif isinstance(wind_direction, np.ndarray):
            wind_direction_array = self._check_wind_direction_array_dims(
                    wind_direction )
            self.set_metadata('WIND_DIRECTION_SOURCE',
                    'numpy.ndarray')

    def get_source_wind(self, eResampleAlg=1, reprojected=True):
        ''' 
        Returns a Nansat object of the auxiliary (model) wind data
        reprojected to the actual SAR image 
        '''

        aux_wind = Nansat(self.get_metadata('WIND_DIRECTION_SOURCE'))
        # Interpolation onto SAR image
        if reprojected:
            aux_wind.reproject(self, eResampleAlg=eResampleAlg, tps=True)

        return aux_wind

    def _get_aux_wind_from_str(self, aux_wind_source):

        try:
            # If a complete filename of wind direction source is given
            aux_wind = Nansat(aux_wind_source)
        except:
            # If only mapper name is given, we add the SAR image
            # timestamp as string. Mappers with this functionality
            # implemented will then find a matching file
            aux_wind = Nansat(aux_wind_source +
                                datetime.strftime(
                                self.SAR_image_time, ':%Y%m%d%H%M'))
        return aux_wind

    def _check_wind_direction_array_dims(self, wind_directions):
        '''
            Check and return the provided array of wind directions
        '''
        if not wind_directions.shape == self.shape():
            raise RuntimeError('The provided wind direction array ' \
                        'must have the same shape as the SAR NRCS')
        return wind_directions

    def _get_wind_direction_array(self, aux_wind):
        '''
            Reproject wind and return the wind directions
        '''
        if not isinstance(aux_wind, Nansat):
            raise ValueError('Wind direction is not available')

        wind_direction_time = aux_wind.get_time()[0]

        # Check time difference between SAR image and wind direction object
        timediff = self.SAR_image_time - wind_direction_time
        try:
            hoursDiff = np.abs(timediff.total_seconds()/3600.)
        except: # for < python2.7
            secondsDiff = (timediff.microseconds +
                            (timediff.seconds + timediff.days *
                            24 * 3600) * 10**6) / 10**6
            hoursDiff = np.abs(secondsDiff/3600.)

        print 'Time difference between SAR image and wind direction: ' \
                + '%.2f' % hoursDiff + ' hours'
        print 'SAR image time: ' + str(self.SAR_image_time)
        print 'Wind dir time: ' + str(wind_direction_time)
        if hoursDiff > 3:
            print '#########################################'
            print 'WARNING: time difference exceeds 3 hours!'
            print '#########################################'

        wind_u_bandNo = aux_wind._get_band_number({
                            'standard_name': 'eastward_wind',
                        })
        wind_v_bandNo = aux_wind._get_band_number({
                            'standard_name': 'northward_wind',
                        })
        # Get wind direction
        u_array = aux_wind[wind_u_bandNo]
        v_array = aux_wind[wind_v_bandNo]
        # 0 degrees meaning wind from North, 90 degrees meaning wind from East
        if u_array==None:
            raise Exception('Could not read wind vectors')
        return np.degrees(np.arctan2(-u_array, -v_array)), \
                wind_direction_time, aux_wind['windspeed']

    def calculate_wind(self, wind_direction='ncep_wind_online',
                        storeModelSpeed=True, eResampleAlg=1):
        '''
            Calculate wind speed from SAR sigma0 in VV polarization
        '''
        wind_direction_array, wind_direction_time, model_windspeed = \
                    self._get_wind_direction_array(
                            self.get_source_wind(eResampleAlg) )

        if isinstance(wind_direction, int) or isinstance(wind_direction, np.ndarray):
            wind_direction_time = None # Not relevant in this case
            storeModelSpeed = False # Not available in this case

        # Calculate SAR wind with CMOD
        # TODO:
        # - add other CMOD versions than CMOD5
        print 'Calculating SAR wind with CMOD...'
        startTime = datetime.now()
        windspeed = cmod5n_inverse(self[self.sigma0_bandNo],
                            np.mod(wind_direction_array -
                                self['SAR_look_direction'], 360),
                            self['incidence_angle'])
        print 'Calculation time: ' + str(datetime.now() - startTime)

        windspeed[np.where(np.isnan(windspeed))] = np.nan
        windspeed[np.where(np.isinf(windspeed))] = np.nan

        ## NOTE: The following will cause a lot of added bands if the function
        #        is repeated - consider changing method to only return data and
        #        make another method that adds the bands

        # Add wind speed and direction as bands
        # TODO: make it possible to update existing bands... See
        # https://github.com/nansencenter/nansat/issues/58
        self.add_band(array=windspeed, parameters={
                        'wkv': 'wind_speed',
                        'name': 'windspeed',
                        'time': self.get_time(self.sigma0_bandNo),
                        'wind_direction_time': wind_direction_time
                })
        self.add_band(array=wind_direction_array, parameters={
                            'wkv': 'wind_from_direction',
                            'name': 'winddirection',
                            'time': wind_direction_time
                })

        if storeModelSpeed:
            self.add_band(array=model_windspeed, parameters={
                            'wkv': 'wind_speed',
                            'name': 'model_windspeed',
                            'time': wind_direction_time,
            })

        # TODO: Replace U and V bands with pixelfunctions
        u = -windspeed*np.sin((180.0 - wind_direction_array)*np.pi/180.0)
        v = windspeed*np.cos((180.0 - wind_direction_array)*np.pi/180.0)
        self.add_band(array=u, parameters={
                            'wkv': 'eastward_wind',
                            'time': wind_direction_time,
        })
        self.add_band(array=v, parameters={
                            'wkv': 'northward_wind',
                            'time': wind_direction_time,
        })

        # set winddir_time to global metadata
        self.set_metadata('winddir_time', str(wind_direction_time))

    def _get_masked_windspeed(self, landmask=True, icemask=True,
            windspeedBand='windspeed'):
        try:
            sar_windspeed = self[windspeedBand]
        except:
            raise ValueError('SAR wind has not been calculated, ' \
                'execute calculate_wind(wind_direction) first.')

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



    def plot(self, filename=None, numVectorsX = 16, show=True,
            clim=[0,20], maskWindAbove=35,
            windspeedBand='windspeed', winddirBand='winddirection',
            northUp_eastRight=True, landmask=True, icemask=True):
        ''' Basic plotting function showing CMOD wind speed
        overlaid vectors in SAR image projection

        parameters
        ----------
        filename : string
        numVectorsX : int
            Number of wind vectors along first dimension
        show : Boolean
        clim : list
            Color limits of the image.
        windspeedBand : string or int
        winddirBand : string or int
        landmask : Boolean
        icemask : Boolean
        maskWindAbove : int

        '''

        try:
            sar_windspeed, palette = self._get_masked_windspeed(landmask,
                    icemask, windspeedBand=windspeedBand)
        except:
            raise ValueError('SAR wind has not been calculated, ' \
                'execute calculate_wind(wind_direction) before plotting.')
        sar_windspeed[sar_windspeed>maskWindAbove] = np.nan

        winddirReductionFactor = np.round(
                self.vrt.dataset.RasterXSize/numVectorsX)

        winddir_relative_up = 360 - self[winddirBand] + \
                                    self.azimuth_y()
        indX = range(0, self.vrt.dataset.RasterXSize, winddirReductionFactor)
        indY = range(0, self.vrt.dataset.RasterYSize, winddirReductionFactor)
        X, Y = np.meshgrid(indX, indY)
        try: # scaling of wind vector length, if model wind is available
            model_windspeed = self['model_windspeed']
            model_windspeed = model_windspeed[Y, X]
        except:
            model_windspeed = 8*np.ones(X.shape)

        Ux = np.sin(np.radians(winddir_relative_up[Y, X]))*model_windspeed
        Vx = np.cos(np.radians(winddir_relative_up[Y, X]))*model_windspeed

        # Make sure North is up, and east is right
        if northUp_eastRight:
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
        plt.clim(clim)
        cbar = plt.colorbar(orientation='horizontal', shrink=.80,
                     aspect=40,
                     fraction=legendFraction, pad=legendPadFraction)
        cbar.ax.set_ylabel('[m/s]', rotation=0) # could replace m/s by units from metadata
        cbar.ax.yaxis.set_label_position('right')
        # TODO: plotting function should be improved to give
        #       nice results for images of all sized
        ax.quiver(X, Y, Ux, Vx, angles='xy', width=0.004,
                    scale=200, scale_units='width',
                    color=[.0, .0, .0], headaxislength=4)
        if filename is not None:
            fig.savefig(filename, pad_inches=0, dpi=dpi)
        if show:
            plt.show()
        return fig

    def save_wind_map_image(self, fileName, scale=None, numArrowsRange=10,
                            landmask=True, colorbar=True, cbar_fontsize=6,
                            drawgrid=True, title=None, title_fontsize=10,
                            edgecolor=None, quiverScaleCriteria=None,
                            tight=True, **kwargs):

        pcolormeshArgs = {'vmin': 0, 'vmax':20}
        for iKey in pcolormeshArgs.keys():
            if iKey in kwargs.keys():
                pcolormeshArgs[iKey] = kwargs.pop(iKey)

        quiverArgs = {'X':None, 'Y':None, 'U':None,
                      'label':None,
                      'labelpos':'E', 'coordinates':'figure',
                      'fontproperties':None, 'width':None}
        popKeys = []
        for iKey in quiverArgs:
            if 'quiver_' + iKey in kwargs.keys():
                quiverArgs[iKey] = kwargs.pop('quiver_'+iKey)
            else:
                popKeys.append(iKey)
        for key in popKeys:
            quiverArgs.pop(key)

        nMap = Nansatmap(self, resolution='l',figsize=(5, 8))

        # use wind direction "to" for calculating u and v
        winddirection = np.mod(self['winddirection'] + 180, 360)
        windspeed = self['windspeed']
        windspeedPcolor = np.copy(windspeed)

        # if data has non-value (dark blue), replace to Nan
        if edgecolor is not None:
            # Replace the edge color (dark blue) to white
            windspeedPcolor[windspeedPcolor == edgecolor] = np.NaN
        # Put colormesh
        nMap.pcolormesh(windspeedPcolor, **pcolormeshArgs)

        # apply landmask to windspeeds
        windspeed = np.ma.masked_where(self.watermask()[1]==2, windspeed)

        # specify the number of quiver
        quiPixelSpacing = int(np.round(windspeed.shape[1]/numArrowsRange))
        numQuiX = int(self.vrt.dataset.RasterXSize / quiPixelSpacing)
        numQuiY = int(self.vrt.dataset.RasterYSize / quiPixelSpacing)

        # compute maximum wind speed on the sea
        maxSpeed = max(windspeed[windspeed.mask == False])
        # compute a scale for quiver lenght
        scale = None
        if quiverScaleCriteria is not None:
            for iKey in quiverScaleCriteria.keys():
                if eval(iKey %maxSpeed):
                    scale = quiverScaleCriteria[iKey]

        # Draw quivers
        Ux = np.sin(np.radians(winddirection)) * windspeed
        Vx = np.cos(np.radians(winddirection)) * windspeed
        nMap.quiver(Ux, Vx, scale=scale,
                    quivectors=(numQuiY, numQuiX), **quiverArgs)

        if colorbar:
            nMap.add_colorbar(fontsize=cbar_fontsize)

        if drawgrid:
            nMap.drawgrid()

        if title is not None:
            plt.title(title, fontsize=title_fontsize)

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
    parser.add_argument('-w', dest='wind_direction',
            default='ncep_wind_online', help='Wind direction filename or constant '
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
        wind_direction = int(args.wind_direction)
    except:
        wind_direction = args.wind_direction

    # Read SAR image
    sw = SARWind(args.SAR_filename, wind_direction, pixelsize=args.pixelsize)

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
