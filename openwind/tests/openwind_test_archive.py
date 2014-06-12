#-------------------------------------------------------------------------------
# Name:         openwind_test_archive.py
# Purpose:      The module contains the paths to files used in the OpenWind test
#               suite
#
# Author:       Morten Wergeland Hansen
# Modified:	Morten Wergeland Hansen
#
# Created:	26.05.2014
# Last modified:12.06.2014 14:07
# Copyright:    (c) NERSC
# License:      
#-------------------------------------------------------------------------------
import os

path = os.path.dirname(__file__)

'''
    Online datasets
'''
asar_agulhas_url = 'ftp://ftp.nersc.no/pub/python_test_data/asar/ASA_WSM_1PNPDE20120327_205532_000002143113_00100_52700_6903.N1'
fname = os.path.basename(asar_agulhas_url)
if not os.path.exists(os.path.join(path,fname)):
    os.system('curl -so ' + os.path.join(path,fname) + ' ' + asar_agulhas_url )
asar_agulhas = os.path.join(path,fname)

ncep_agulhas_url = 'ftp://ftp.nersc.no/pub/python_test_data/ncep/gfs/gfs20120328/gfs.t00z.master.grbf00'
fname = os.path.basename(ncep_agulhas_url)
if not os.path.exists(os.path.join(path,fname)):
    os.system('curl -so ' + os.path.join(path,fname) + ' ' + ncep_agulhas_url )
ncep_agulhas = os.path.join(path,fname)


'''
    Local datasets (e.g. Radarsat-2)
'''
try:
    from openwind_local_archive import *
except ImportError as e:
    print 'ImportError: {0}'.format(e)
    print '\nPlease place a module openwind_local_archive on your python ' \
        'path and add paths to local test data if you have any, then ' \
        'update openwind.tests.test_sarwind.py.\n'
except IOError as e:
    print e.message



'''
    Check existence
'''
if not os.path.isfile(asar_agulhas) or not os.path.isfile(ncep_agulhas):
    asar_agulhas = None
    ncep_agulhas = None

'''
    Test data should be contained in an instance of the TestData class so we're
    able to correctly remove downloaded files after the tests
'''
class TestData():
    asar = []
    radarsat2 = []
    ncep4asar = []
    noData = True

    def __init__(self):
        # OBS: SAR and wind data must be added in pairs for each test
        if asar_agulhas:
            self.asar.append(asar_agulhas)
        if ncep_agulhas:
            self.ncep4asar.append(ncep_agulhas)
        if 'rs2' in globals():
            self.radarsat2.append(rs2)
        if self.asar:
            self.noData = False
        if self.radarsat2:
            self.noData = False

    def __del__(self):
        '''
            Delete any downloaded files
        '''
        for a in self.asar:
            if os.path.isfile(a):
                os.unlink(a)
        for n in self.ncep4asar:
            if os.path.isfile(n):
                os.unlink(n)
        self.asar = []
        self.ncep4asar = []
        self.noData = True
