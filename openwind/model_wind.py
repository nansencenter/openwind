#-------------------------------------------------------------------------------
# Name:		model_wind.py
# Purpose:      
#
# Author:       Morten Wergeland Hansen
# Modified:	Morten Wergeland Hansen
#
# Created:	15.08.2013
# Last modified:24.10.2013 14:29
# Copyright:    (c) NERSC
# License:      GNU GPL
#-------------------------------------------------------------------------------
from nansat import Nansat

import os
import numpy as np

import pdb

# Locations of the grib-files are hardcoded to my local directory for now.
# Later, the model wind field should be retrieved from an online service...

# REMEMBER TO MODIFY THESE LISTS TO SUIT YOUR LOCAL SYSTEM:
HIRLAM_wind_dir = [ # Add the location of your hirlam dir
                    '/Volumes/sat/auxdata/model/met.no/',
                    '/Users/mortenh/met.no/',
                ]
NCEP_wind_dir = [   # Add the location of your ncep dir
                    '/Volumes/sat/auxdata/model/ncep/',
                    '/Users/mortenh/ncep/',
                ]

class ModelWind(Nansat, object):

    def __init__(self, sar, *args, **kwargs):
        '''
            Set parameters needed to get the model wind field

            Presently, a Nansat sar object is used. This is a bit heavy and
            should probably be replaced by a simpler nansat.Domain object and
            a time.
        '''
        if isinstance(sar,basestring):
            self.sar = Nansat(sar)
        else:
            self.sar = sar
        gribfile = self.get_hirlam_wind()
        success = True
        if os.path.exists(gribfile):
            super(ModelWind, self).__init__(gribfile, *args, **kwargs)
            if not self.check_coverage():
                ## delete HIRLAM file
                #os.remove(self.windGribFile)
                #super(ModelWind, self).__del__()
                print 'No HIRLAM coverage, trying NCEP...'
                gribfile = self.get_ncep_wind()
                success = False
        else:
            print 'No HIRLAM wind file, trying NCEP...'
            gribfile = self.get_ncep_wind()
            success = False
        if not success:
            if os.path.exists(gribfile):
                success = True
                super(ModelWind, self).__init__(gribfile, *args, **kwargs)
            else:
                raise IOError('Cannot locate model wind field information')
        # Reproject model wind field to SAR image coverage
        self.reproject(self.sar)

    def __del__(self):
        if hasattr(self, 'windGribFile') and os.path.exists(self.windGribFile):
            # Delete grib-file
            os.remove(self.windGribFile)

    def check_coverage(self):
        lon_covered =\
                np.min(self.sar.get_corners()[0])>=np.min(self.get_corners()[0]) \
                and \
                np.max(self.sar.get_corners()[0])<=np.max(self.get_corners()[0])
        lat_covered =\
                np.min(self.sar.get_corners()[1])>=np.min(self.get_corners()[1]) \
                and \
                np.max(self.sar.get_corners()[1])<=np.max(self.get_corners()[1])
        if not lon_covered or not lat_covered:
            return False
        else:
            return True

    def get_hirlam_wind(self):
        #'windfile': 'model_wind/HIRLAM_10kmEurope_20120615_1200.grib'
        # HIRLAM wind file 
        estr = 'Make sure one of these directories are connected:\n\n'
        for dir in HIRLAM_wind_dir:
            estr = estr + dir + '\n'
            if os.path.exists(dir):
                hdir = dir
                break

        if not 'hdir' in locals():
            raise Exception, estr

        hour = '1200'
        if self.sar.get_time()[0].hour + self.sar.get_time()[0].minute/60.0 < 13.5:
            hour = '0000'
        hfile = os.path.join( hdir, 'HIRLAM_10kmEurope_' + 
                    self.sar.get_time()[0].strftime('%Y%m%d') + '_' + hour + '.grib' )

        print hfile
        return hfile

    def get_ncep_wind(self):
        # NCEP wind file 
        estr = 'Make sure one of these directories are connected:\n\n'
        for dir in NCEP_wind_dir:
            estr = estr + dir + '\n'
            if os.path.exists(dir):
                ndir = dir
                break
        if not 'ndir' in locals():
            raise Exception, estr

        basehour = np.floor((self.sar.get_time()[0].hour + self.sar.get_time()[0].minute/60.0 +
            3.0/2.0 )/6.0)*6.0
        basehour = np.min([18, basehour])
        if self.sar.get_time()[0].hour + self.sar.get_time()[0].minute/60.0 - basehour > 1.5:
            forecasthour = 3
        else:
            forecasthour = 0
        nfile = os.path.join( ndir, 'gfs',
                'gfs'+self.sar.get_time()[0].strftime('%Y%m%d'),
                'gfs.t%.2dz.master.grbf%.2d' %(basehour, forecasthour) )
        print nfile
        return nfile
