# Name:		open_wind.py
# Purpose:      
# Authors:      Morten Wergeland Hansen (mortenwh)
# License:      This file is part of OPENWIND. You can redistribute it or
#               modify under the terms of GNU General Public License, v.3
#               http://www.gnu.org/licenses/gpl-3.0.html
from nansat import Nansat
import numpy as np

from model_wind import ModelWind
from cmod5n import cmod5n_inverse

import scipy

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
        super(SARWind, self).__init__(*args,**kwargs)

        # Check that this is a SAR acquisition that has a NanSat mapper and
        # that it has VV pol nrcs
        names = [self.bands()[k]['name'] for k in self.bands().keys()] 
        if not 'sigma0_VV' in names:
            raise TypeError(self.filename + ' is not a valid NanSat SAR image file')

        # Get model wind field
        model_wind = ModelWind(self)
        winddir = model_wind.get_GDALRasterBand('winddirection').ReadAsArray()

        # Note that this is the look direction in the center of the domain. For
        # longer domains, especially at high latitudes, the look direction
        # may vary a lot over the domain, and using the center angle will be a
        # coarse approximation.
        look_direction = float(self.vrt.dataset.GetMetadataItem('SAR_center_look_direction'))
        winddir_relative = np.mod(winddir - look_direction, 360)
        windspeed = self.get_cmod_wind(winddir_relative)

        u = windspeed*np.sin((winddir+180)/np.pi)
        v = windspeed*np.cos((winddir+180)/np.pi)

        self.add_band(array=u, parameters={
                        'wkv': 'eastward_wind',
                    })
        self.add_band(array=v, parameters={
                        'wkv': 'northward_wind',
                    })

        names = [self.bands()[k]['name'] for k in self.bands().keys()] 
        # Add pixelfunctions for retrieving wind speed and direction
        metaDict = [
                    {
                        'src': [{
                                'SourceFilename': self.raw.fileName,
                                'SourceBand': 
                                    int(np.where(np.array(names)=='U')[0][0])
                            },{
                                'SourceFilename': self.raw.fileName,
                                'SourceBand':
                                    int(np.where(np.array(names)=='V')[0][0])
                                }],
                        'dst': {
                                'wkv': 'wind_from_direction',
                                'name': 'winddirection',
                                'PixelFunctionType': 'UVToDirectionFrom',
                            }
                    },{
                        'src': [{
                                'SourceFilename': self.raw.fileName,
                                'SourceBand': 
                                    int(np.where(np.array(names)=='U')[0][0])
                            },{
                                'SourceFilename': self.raw.fileName,
                                'SourceBand':
                                    int(np.where(np.array(names)=='V')[0][0])
                                }],
                        'dst': {
                            'wkv': 'wind_speed',
                            'name': 'windspeed',
                            'PixelFunctionType': 'UVToMagnitude',
                            }
                    },
                ]
        self.raw._create_bands(metaDict)
        # Add metadata:
        #               - Model wind field
        #               - CMOD-function
        #               - etc
        self.raw.dataset.SetMetadataItem('Model wind field',
                model_wind.mapper[7:])

        # copy raw VRT object to the current vrt
        self.vrt = self.raw.copy()

    def invalid2nan(self, bandName):
        nparr = self.get_GDALRasterBand(bandName).ReadAsArray()
        if not np.isrealobj(nparr[0][0]):
            nparr = np.abs(nparr)
        if self.get_metadata(bandID=bandName).has_key('_FillValue'):
            fillValue = float(self.get_metadata(bandID=bandName)['_FillValue'])
            # Set invalid and missing data to np.NaN
            # This is relevant for nansat data retrieved from netcdf
            nparr[np.where(nparr==fillValue)]=np.nan
        return nparr

    def get_cmod_wind(self, winddir_relative):
        s0 = self.invalid2nan('sigma0_VV')
        # consider changing cmod to use only one line of incidence angles
        windspeedcmod5 = cmod5n_inverse( s0, winddir_relative,
                self.incidence_angle() )
        windspeedcmod5[np.where(np.isnan(s0))] = np.nan
        windspeedcmod5[np.where(np.isinf(s0))] = np.nan
        return windspeedcmod5

    def incidence_angle(self, **kwargs):
        inci = self.get_GDALRasterBand('incidence_angle').ReadAsArray(**kwargs)
        if not np.isrealobj(inci[0][0]):
            # try to open non-complex band
            if self.has_band('incidence_angle_noncomplex'):
                inci = self.nObj.get_GDALRasterBand('incidence_angle_noncomplex').ReadAsArray(**kwargs)
            else:
                inci = np.abs(self.nObj.get_GDALRasterBand('incidence_angle').ReadAsArray(**kwargs))
                self.nObj.add_band(array=inci, parameters={
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

    def plot_example(self):
        resize_factor = 0.2
        resampleAlg = 0 # nearest neighbour
        look_direction = float(self.get_metadata('SAR_center_look_direction'))
        self.resize(resize_factor, eResampleAlg=resampleAlg)
        speed = self.invalid2nan('windspeed')
        dirGeo = self.get_GDALRasterBand('winddirection').ReadAsArray()
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
        Q = map.barbs(x[dd-1::dd,::dd], y[dd-1::dd,::dd],
                self.get_GDALRasterBand('U').ReadAsArray()[dd-1::dd,::dd],
                self.get_GDALRasterBand('V').ReadAsArray()[dd-1::dd,::dd])
        #Quiver
        #Q = map.quiver(x[::dd,dd-1::dd], y[::dd,dd-1::dd],
        #        self.get_GDALRasterBand('east_wind').ReadAsArray()[::dd,dd-1::dd],
        #        self.get_GDALRasterBand('north_wind').ReadAsArray()[::dd,dd-1::dd],
        #        width=0.001)
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
        fig.savefig( 'windfield'+os.path.basename(self.fileName)+'.png', facecolor='w', edgecolor='w', dpi=300,
                bbox_inches="tight", pad_inches=0.1)

        self.resize()


