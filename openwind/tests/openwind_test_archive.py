#-------------------------------------------------------------------------------
# Name:         openwind_test_archive.py
# Purpose:      The module contains the paths to files used in the OpenWind test
#               suite
#
# Author:       Morten Wergeland Hansen
# Modified:	Morten Wergeland Hansen
#
# Created:	26.05.2014
# Last modified:07.07.2014 15:31
# Copyright:    (c) NERSC
# License:      
#-------------------------------------------------------------------------------
import os

from nansat.tests.nansat_test_archive import TestData
dirname_test_data = os.path.join(
                        os.path.dirname(os.path.abspath(__file__)),
                        'test_data')
if not os.path.exists(dirname_test_data):
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
class OpenWindTestData(TestData):

    radarsat2 = []
    ncep4asar = []

    def __init__(self):
        # OBS: SAR and wind data must be added in parallel
        super(OpenWindTestData, self).get_asar_agulhas()
        self.get_ncep_agulhas()
        self.get_local()
        if self.asar:
            self.noData = False
        if self.radarsat2:
            self.noData = False

    def get_ncep_agulhas(self):
        ncep_agulhas_url = 'ftp://ftp.nersc.no/pub/python_test_data/ncep/gfs/gfs20120328/gfs.t00z.master.grbf00'
        fname = os.path.basename(ncep_agulhas_url)
        ncep_agulhas = os.path.join(dirname_test_data,fname)
        if not os.path.exists(ncep_agulhas):
            os.system('curl -so ' + ncep_agulhas + ' ' + ncep_agulhas_url )
        if not os.path.isfile(ncep_agulhas):
            ncep_agulhas = None
            warnings.warn( "Could not access ftp-site with test data - contact " \
                    "morten.stette@nersc.no to get the ftp-server at NERSC restarted" )
        else:
            self.ncep4asar.append(ncep_agulhas)

    def get_local(self):
        # check for data from local archives
        if 'rs2' in globals() and rs2:
            self.radarsat2.append(rs2)
        if 'rs2_quad' in globals() and rs2_quad:
            self.radarsat2.append(rs2_quad)

