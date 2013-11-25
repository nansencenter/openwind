#!/usr/bin/env python
#
# Utility to download GRIB files with 
# NCEP GFS model wind from NOMADS web servers


import sys
import os
from datetime import datetime, timedelta
from nansat import Nansat

openWindFolder = os.path.dirname(__file__) + '/'
get_inv = openWindFolder + 'get_inv.pl '
get_grib = openWindFolder + 'get_grib.pl '


def download_ncep(time, outFolder = ''):

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


############################
# Command line usage
############################
if __name__ == "__main__":
        
    try:
        time = datetime.strptime(sys.argv[1], '%Y%m%d%H%M')
    except:
        print 'Usage: ' 
        print sys.argv[0] + ' <time>'
        print sys.argv[0] + ' %Y%m%d%H%M'
        print sys.argv[0] + ' 201311070600'
        print ''
        sys.exit()

    fileName = download_ncep(time)
    if fileName is not None:
        print 'Downloaded file: ' + fileName
    else:
        print 'NCEP field not available at ' + str(time)
