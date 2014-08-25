#!/usr/bin/env python
#-------------------------------------------------------------------------------
# Name:		test_sarwind.py
# Purpose:
#
# Author:       Morten Wergeland Hansen
# Modified:	Morten Wergeland Hansen
#
# Created:	26.05.2014
# Last modified:14.07.2014 13:27
# Copyright:    (c) NERSC
# License:
#-------------------------------------------------------------------------------
import unittest
import filecmp
import os, sys
import numpy as np

from openwind import SARWind
from nansat.nansat import Nansat, Domain

import openwind.tests.openwind_test_archive as ota

dirname_test_plots = os.path.join(
                        os.path.dirname(os.path.abspath(__file__)), 'plots')

'''
Downloads online test data, saves it locally and performs tests.

This test suite does not test command line usage. See
http://hg.python.org/cpython/file/default/Lib/test/test_cmd_line_script.py for
several examples of how to test handling of command line arguments.
'''

class SARWindTest(unittest.TestCase):
    def setUp(self):
        self.test_data = ota.OpenWindTestData()
        if self.test_data.noData:
            raise ValueError('No test data available')

    def test_sarwind_using_default(self):
        if len(self.test_data.radarsat2)==0:
            raise IOError('No Radarsat-2 data - try adding some as ' \
                    'described in templates/openwind_local_archive.py' )
        for i in range(len(self.test_data.radarsat2)):
            w = SARWind(self.test_data.radarsat2[i])
        if sys.version_info < (2, 7):
            type(w) == SARWind
        else:
            self.assertIsInstance(w, SARWind)

    def test_sarwind_using_filenames(self):
        if len(self.test_data.asar)==0:
            raise IOError('No ASAR data - try adding some as ' \
                    'described in templates/openwind_local_archive.py' )
        for i in range(len(self.test_data.asar)):
            w = SARWind(self.test_data.asar[i],
                    wind_direction=self.test_data.ncep4asar[i])
        if sys.version_info < (2, 7):
            type(w) == SARWind
        else:
            self.assertIsInstance(w, SARWind)

    def test_sarwind_using_asar_filename_ncep_nansat(self):
        if len(self.test_data.asar)==0:
            raise IOError('No ASAR data - try adding some as ' \
                    'described in templates/openwind_local_archive.py' )
        for i in range(len(self.test_data.asar)):
            mw = Nansat(self.test_data.ncep4asar[i])
            w = SARWind(self.test_data.asar[i], wind_direction=mw)
        if sys.version_info < (2, 7):
            type(w) == SARWind
        else:
            self.assertIsInstance(w, SARWind)

    def test_sarwind_using_asar_nansat_ncep_filename(self):
        if len(self.test_data.asar)==0:
            raise IOError('No ASAR data - try adding some as ' \
                    'described in templates/openwind_local_archive.py' )
        for i in range(len(self.test_data.asar)):
            asar = Nansat(self.test_data.asar[i])
            w = SARWind(asar, wind_direction=self.test_data.ncep4asar[i])
        if sys.version_info < (2, 7):
            type(w) == SARWind
        else:
            self.assertIsInstance(w, SARWind)

    def test_sarwind_using_asar_nansat_ncep_nansat(self):
        if len(self.test_data.asar)==0:
            raise IOError('No ASAR data - try adding some as ' \
                    'described in templates/openwind_local_archive.py' )
        for i in range(len(self.test_data.asar)):
            asar = Nansat(self.test_data.asar[i])
            mw = Nansat(self.test_data.ncep4asar[i])
            w = SARWind(asar, wind_direction=mw)
        if sys.version_info < (2, 7):
            type(w) == SARWind
        else:
            self.assertIsInstance(w, SARWind)

    def test_nansat_reproject(self):
        if len(self.test_data.asar)==0:
            raise IOError('No ASAR data - try adding some as ' \
                    'described in templates/openwind_local_archive.py' )
        asar = Nansat(self.test_data.asar[0])
        asar.resize(pixelsize=500, eResampleAlg=1)
        mw = Nansat(self.test_data.ncep4asar[0])
        mw.reproject(asar, eResampleAlg=1)
        if sys.version_info < (2, 7):
            type(mw[1]) == np.ndarray
        else:
            self.assertIsInstance(mw[1], np.ndarray)

    def test_plot_agulhas(self):
        '''
        Tests the SARWind.plot method. Remember to update the file
        basic_agulhas_ref.png if plot method is changed.

        If adding more test plots, please make sure the indices are correct
        '''
        if len(self.test_data.asar)==0:
            raise IOError('No ASAR data - try adding some as ' \
                    'described in templates/openwind_local_archive.py' )
        w = SARWind(self.test_data.asar[0],
                wind_direction=self.test_data.ncep4asar[0])
        w.plot(filename=os.path.join(dirname_test_plots,
            'agulhas_test_plot.png'), show=False, landmask=False)

        if sys.version_info < (2, 7):
            self.assertTrue(filecmp.cmp(
                os.path.join(dirname_test_plots,'agulhas_test_plot.png'),
                os.path.join(dirname_test_plots,'agulhas_test_plot_ref.png'))
                in [True])
        else:
            self.assertTrue(filecmp.cmp(
                os.path.join(dirname_test_plots,'agulhas_test_plot.png'),
                os.path.join(dirname_test_plots,'agulhas_test_plot_ref.png')))

    def test_plot_barents(self):
        if len(self.test_data.radarsat2)==0:
            raise IOError('No Radarsat-2 data - try adding some as ' \
                    'described in templates/openwind_local_archive.py' )
        w = SARWind(self.test_data.radarsat2[0])
        w.plot(filename=os.path.join(dirname_test_plots,
            'rs2_barents_test_plot.png'), show=False, landmask=False)
        if sys.version_info < (2, 7):
            self.assertTrue(filecmp.cmp(
                os.path.join(dirname_test_plots, 'rs2_barents_test_plot.png'),
                os.path.join(dirname_test_plots, 'rs2_barents_test_plot_ref.png'))
                in [True])
        else:
            self.assertTrue(filecmp.cmp(
                os.path.join(dirname_test_plots,'rs2_barents_test_plot.png'),
                os.path.join(dirname_test_plots,'rs2_barents_test_plot_ref.png')))

    def tearDown(self):
        # Delete test plots
        if os.path.exists( os.path.join( dirname_test_plots,
                'agulhas_test_plot.png' ) ):
            os.unlink(os.path.join(dirname_test_plots,'agulhas_test_plot.png'))

        if os.path.exists( os.path.join( dirname_test_plots,
                'rs2_barents_test_plot.png' ) ):
            os.unlink(os.path.join(dirname_test_plots,'rs2_barents_test_plot.png'))


if __name__=='__main__':
    unittest.main()
