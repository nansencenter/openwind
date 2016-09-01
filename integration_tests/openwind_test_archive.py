#-------------------------------------------------------------------------------
# Name:         openwind_test_archive.py
# Purpose:      The module contains the paths to files used in the OpenWind test
#               suite
#
# Author:       Morten Wergeland Hansen
# Modified:	Morten Wergeland Hansen
#
# Created:	26.05.2014
# Last modified:21.11.2014 10:33
# Copyright:    (c) NERSC
# License:      
#-------------------------------------------------------------------------------
import os
import warnings

dirname_test_data = os.path.join(os.path.expanduser('~'), 'openwind_test_data')
if not os.path.isdir(dirname_test_data):
    os.mkdir(dirname_test_data)

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
    Test data should be contained in an instance of the TestData class so we're
    able to correctly remove downloaded files after the tests
'''
class OpenWindTestData():

    asar = {}
    radarsat2 = {}
    sentinel1a = {}
    ncep4asar = {}
    noData = True

    def __init__(self):
        # OBS: SAR and wind data must be added in parallel
        self.get_asar_agulhas()
        self.get_ncep_agulhas()
        self.get_local()
        if self.asar:
            self.noData = False
        if self.radarsat2:
            self.noData = False

    def get_ncep_agulhas(self):
        ncep_agulhas_url = 'ftp://ftp.nersc.no/pub/nansat/test_data/ncep/gfs20120328.t00z.master.grbf00'
        fname = os.path.basename(ncep_agulhas_url)
        ncep_agulhas = os.path.join(dirname_test_data,fname)
        if not os.path.isfile(ncep_agulhas):
            os.system('curl -so ' + ncep_agulhas + ' ' + ncep_agulhas_url )
        if not os.path.isfile(ncep_agulhas):
            ncep_agulhas = None
            warnings.warn( "Could not access ftp-site with test data - contact " \
                    "morten.stette@nersc.no to get the ftp-server at NERSC restarted" )
        else:
            self.ncep4asar['agulhas'] = ncep_agulhas

    def get_local(self):
        # check for data from local archives
        if 'rs2' in globals() and rs2:
            self.radarsat2.append(rs2)
        if 'rs2_quad' in globals() and rs2_quad:
            self.radarsat2.append(rs2_quad)

    def get_asar_agulhas(self):
        asar_agulhas_url = 'ftp://ftp.nersc.no/pub/nansat/test_data/asar/ASA_WSM_1PNPDE20120327_205532_000002143113_00100_52700_6903.N1'
        fname = os.path.basename(asar_agulhas_url)
        asar_agulhas = os.path.join(dirname_test_data, fname)
        if not os.path.isfile(asar_agulhas):
            os.system('curl -so ' + asar_agulhas + ' ' + asar_agulhas_url )
        if not os.path.isfile(asar_agulhas):
            asar_agulhas = None
            warnings.warn( "Could not access ftp-site with test data - contact " \
                    "morten.stette@nersc.no to get the ftp-server at NERSC restarted" )
        else:
            self.asar['agulhas'] = asar_agulhas

