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


import matplotlib 
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import geographiclib.geodesic
import mpl_toolkits
from mpl_toolkits.basemap import Basemap
import figure # figure module in nansat

import pdb

class SARWind(Nansat, object):
    '''
    A class for retrieving wind speed from SAR images using cmod
    '''

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
            raise TypeError(self.fileName + ' is not a valid NanSat SAR image file')
        
        if os.path.splitext(self.fileName)[1]=='.nc':
            return

        # Reduce size of image to 500 m pixels in range
        line_spacing = float(self.get_metadata()['LINE_SPACING'])
        self.resize(line_spacing/500.)
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

        self.add_band(array=windspeed, parameters={
                        'wkv': 'wind_speed',
                        'name': 'windspeed',
                    })
        self.add_band(array=winddir, parameters={
                        'wkv': 'wind_from_direction',
                        'name': 'winddirection',
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
        s0 = np.abs(self['sigma0_VV'])
        # consider changing cmod to use only one line of incidence angles
        windspeedcmod5 = cmod5n_inverse( s0, winddir_relative,
                np.abs(self.incidence_angle()) )
        windspeedcmod5[np.where(np.isnan(s0))] = np.nan
        windspeedcmod5[np.where(np.isinf(s0))] = np.nan
        return windspeedcmod5

    def incidence_angle(self, **kwargs):
        inci = self.get_GDALRasterBand('incidence_angle').ReadAsArray(**kwargs)
        if not np.isrealobj(inci[0][0]):
            # try to open non-complex band
            if self.has_band('incidence_angle_noncomplex'):
                inci = self.get_GDALRasterBand('incidence_angle_noncomplex').ReadAsArray(**kwargs)
                # Set invalid and missing data to np.nan
                if inci.GetMetadata().has_key('_FillValue'):
                    fillValue = float(inci.GetMetadata()['_FillValue'])
                    inci[np.where(inci==fillValue)]=np.nan
                inci[np.where(np.isinf(inci))] = np.nan
            else:
                inci = np.abs(self.get_GDALRasterBand('incidence_angle').ReadAsArray(**kwargs))
                self.add_band(array=inci, parameters={
                        'name': 'incidence_angle_noncomplex',
                        'long_name': 'Non-complex incidence angle for correct netcdf export',
                        'wkv': 'angle_of_incidence',
                        'dataType': 6,
                    }) 
        ind=np.where(inci==0)
        # set incidence angles which are 0 to np.nan (could also make some
        # interpolation but that is for later improvements...)
        inci[ind] = np.nan
        return inci

    def plot_example(self, model_wind=None, filename='windfield.png', dir='.'):

        uu = self['U']
        vv = self['V']
        if model_wind:
            uu = model_wind['U']
            vv = model_wind['V']
        look_direction = float(self.get_metadata('SAR_center_look_direction'))
        speed = self['windspeed']
        dirGeo = self['winddirection']
        dirLookRelative = np.mod(np.subtract( dirGeo, look_direction ), 360)
        dirRange = -np.sin(dirLookRelative*np.pi/180.)
        dirAzim = np.cos(dirLookRelative*np.pi/180.)
        x=np.arange(np.shape(dirRange)[0])
        y=np.arange(np.shape(dirRange)[1])
        X, Y = np.meshgrid(y, x)
        wsmin = scipy.stats.nanmedian(speed, axis=None) - 2*scipy.stats.nanstd(speed,axis=None)
        if wsmin<0:
            wsmin = 0
        wsmax = scipy.stats.nanmedian(speed, axis=None) + 2*scipy.stats.nanstd(speed,axis=None)
        #wsmin = 2
        #wsmax = 20

        lonMax = np.max(self.get_corners()[0])
        lonMean = np.mean(self.get_corners()[0])
        lonMin = np.min(self.get_corners()[0])
        latMax = np.max(self.get_corners()[1])
        latMean = np.mean(self.get_corners()[1])
        latMin = np.min(self.get_corners()[1])
        width = geographiclib.geodesic.Geodesic.WGS84.Inverse(latMean,
                lonMin, latMean, lonMax)
        height = geographiclib.geodesic.Geodesic.WGS84.Inverse(latMin,
                lonMean, latMax, lonMean)

        #fig, ax = self.mkfig()
        # make figsize dynamic
        dpi = 300.
        ysize = self.shape()[0]/dpi
        xsize = self.shape()[1]/dpi
        #print 'xsize: %d\nysize: %d' %(xsize,ysize)
        titleSize = 1.5*xsize #self.shape()[1]/500
        legendSize = xsize #*2.5/3
        tickSize = xsize*2/3 #self.shape()[1]/1000
        fig = plt.figure(num=1, figsize=(xsize,ysize), dpi=dpi)
        plt.clf()
        map = Basemap(projection='laea', 
                        width=width['s12'],
                        height=height['s12'],
                        lat_ts = np.mean(self.get_corners()[1]),
                        lon_0 = np.mean(self.get_corners()[0]),
                        lat_0 = np.mean(self.get_corners()[1]),
                        resolution='i')

        lon,lat=self.get_geolocation_grids()
        x,y = map(lon,lat)
        mappable = map.pcolormesh(x,y,speed,vmin=wsmin,vmax=wsmax,shading='flat',cmap=plt.cm.jet)

        dd = int(np.round(np.shape(speed)[1]/10))
        # Meteorological barbs
        Q = map.barbs(x[dd::dd-1,::dd], y[dd::dd-1,::dd], uu[dd::dd-1,::dd],
                vv[dd::dd-1,::dd])
        map.fillcontinents(color='#cc9966',lake_color='#99ffff')

        nLines = 3.
        map.drawparallels(np.arange(lat.min()+(lat.max()-lat.min())/(nLines*2),lat.max(),
            (lat.max()-lat.min())/nLines), 
            linewidth=0.2,
            labels=[True,False,False,False], fontsize=tickSize, fmt='%.1f')
        map.drawmeridians(np.arange(lon.min()+(lon.max()-lon.min())/(nLines*2),lon.max(),
            (lon.max()-lon.min())/nLines),
            linewidth=0.2,
            labels=[False,False,False,True], fontsize=tickSize, fmt='%.1f')

        #cb = self.add_cbar(fig,ax, mappable)
        cb = fig.colorbar(mappable,orientation='horizontal', pad=0.05)
        cb.set_label('m/s',fontsize=legendSize)

        ax = plt.gca()
        plt.axes(cb.ax)
        plt.xticks(fontsize=tickSize)
        plt.axes(ax)

        #if wtype=='model':
        #    t=plt.title('Model wind field',fontsize=titleSize) # Doesn't work: ,pad_inches=0.1)
        #elif wtype=='cmod':
        #    t=plt.title('CMOD wind field',fontsize=titleSize)
        #t.set_position((0.5,1.02))
        fig.savefig( os.path.join(dir,filename), facecolor='w', edgecolor='w', dpi=300,
                bbox_inches="tight", pad_inches=0.1)
