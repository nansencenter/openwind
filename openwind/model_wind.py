#-------------------------------------------------------------------------------
# Name:		model_wind.py
# Purpose:      
#
# Author:       Morten Wergeland Hansen, Knut-Frode Dagestad
# Modified:	Morten Wergeland Hansen
#
# Created:	15.08.2013
# Last modified:19.12.2013 16:43
# Copyright:    (c) NERSC
# License:      GNU GPL
#-------------------------------------------------------------------------------
from nansat import Nansat

import os
import numpy as np
from datetime import datetime, timedelta

import pdb

openWindFolder = os.path.dirname(__file__)
get_inv = os.path.join(openWindFolder, 'get_inv.pl ')
get_grib = os.path.join(openWindFolder, 'get_grib.pl ')

class ModelWind(Nansat, object):

    def __init__(self, wind=None, time=None, domain=None, *args, **kwargs):
        '''
            Get model wind field, and optionally reproject it to a given
            domain.

            Either filename or Nansat object, or time and domain is required
            as input.
            
            Parameters
            -----------
            wind : string or Nansat object
            time : datetime
                the appreciated time of the model wind field
            domain : nansat Domain
        '''
        if not time and not wind:
            raise ValueError(
                    "At least 'time' or 'wind' is required as input"
                    )
        
        self.downloaded = ''

        if wind:
            if isinstance(wind, str):
                super(ModelWind, self).__init__(wind)
            elif isinstance(wind, Nansat):
                super(ModelWind, self).__init__(wind.fileName)
                self.reproject(wind)
            # Check if Nansat object contains wind direction
            try:
                wind_u_bandNo = self._get_band_number(
                                {'standard_name': 'eastward_wind'})
            except:
                raise TypeError(self.fileName +
                    ' does not contain wind direction')
        else:
            if not isinstance(time, datetime):
                raise TypeError('Time input must be a datetime object')
            else:
                self.downloaded = self.download_ncep(time)
                if self.downloaded is None:
                    raise ValueError('No NCEP file available for time of ' \
                        'SAR image. Can not calculate SAR wind speed.')
                super(ModelWind, self).__init__(self.downloaded)
                wind_u_bandNo = self._get_band_number(
                                {'standard_name': 'eastward_wind'})

        self.time = self.get_time(wind_u_bandNo)

        if domain:
            # Bi-linear interpolation onto given domain
            self.reproject(domain, eResampleAlg=1)

    def __del__(self):
        if self.downloaded and os.path.exists(self.downloaded):
            # Delete downloaded wind file
            os.remove(self.downloaded)

    def check_coverage(self, domain):
        lon_covered =\
                np.min(domain.get_corners()[0])>=np.min(self.get_corners()[0]) \
                and \
                np.max(domain.get_corners()[0])<=np.max(self.get_corners()[0])
        lat_covered =\
                np.min(domain.get_corners()[1])>=np.min(self.get_corners()[1]) \
                and \
                np.max(domain.get_corners()[1])<=np.max(self.get_corners()[1])
        if not lon_covered or not lat_covered:
            return False
        else:
            return True

    def download_ncep(self, time, outFolder = ''):
    
        time = time.replace(tzinfo=None) # Remove timezone information
        # Find closest 6 hourly modelrun and forecast hour
        modelRunHour = round((time.hour + time.minute/60.)/6)*6
        nearestModelRun = datetime(time.year, time.month, time.day) \
            + timedelta(hours=modelRunHour)
        forecastHour = (time - nearestModelRun).total_seconds()/3600.
        if forecastHour < 1.5:
            forecastHour = 0
        else:
            forecastHour = 3
    
        # Try to get NRT data from 
        # ftp://ftp.ncep.noaa.gov/pub/data/nccf/com/gfs/prod/
        # Avaliable approximately the latest month
        url = 'ftp://ftp.ncep.noaa.gov/pub/data/nccf/com/gfs/prod/' \
                + 'gfs.' + time.strftime('%Y%m%d') + '%.2d' % modelRunHour \
                + '/gfs.t' + '%.2d' % modelRunHour + 'z.master.grbf' \
                + '%.2d' % forecastHour + '.10m.uv.grib2'
        outFileName = outFolder + 'ncep_gfs_' + nearestModelRun.strftime(
                    '%Y%m%d_%HH_') + '%.2d' % forecastHour + '.10m.uv.grib2'
        if os.path.exists(outFileName):
            print 'NCEP wind is already downloaded: ' + outFileName
            return outFileName
        else:
            os.system('curl -o ' + outFileName + ' ' + url)
            if os.path.exists(outFileName):
                print 'Downloaded ' + outFileName
                return outFileName
            else:
                print 'NRT GRIB file not available: ' + url
    
        # If NRT file not available, 
        # coninue to search for archived file
        url = 'http://nomads.ncdc.noaa.gov/data/gfs4/' + \
            nearestModelRun.strftime('%Y%m/%Y%m%d/')
        baseName = 'gfs_4_' + nearestModelRun.strftime('%Y%m%d_') \
                + nearestModelRun.strftime('%H%M_') \
                + '%.3d' % forecastHour
        fileName = baseName + '.grb2'
        outFileName = outFolder + fileName
        print 'Downloading ' + url + fileName
    
        # Download subset of grib file
        if not os.path.exists(fileName):
            command = get_inv + url + baseName \
                  + '.inv | egrep "(:UGRD:10 m |:VGRD:10 m )" | ' \
                  + get_grib + url + fileName + ' ' + outFileName
            os.system(command)
        else:
            print 'Already downloaded'
        if os.path.exists(outFileName):
            return outFileName
        else:
            return None

