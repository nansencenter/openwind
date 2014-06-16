#!/usr/bin/env python
#-------------------------------------------------------------------------------
# Name:		test_sarwind.py
# Purpose:      
#
# Author:       Morten Wergeland Hansen
# Modified:	Morten Wergeland Hansen
#
# Created:	26.05.2014
# Last modified:16.06.2014 11:54
# Copyright:    (c) NERSC
# License:      
#-------------------------------------------------------------------------------
import unittest
import filecmp
import os
dirname = os.path.dirname(__file__)

from openwind import *
from nansat import *

import openwind_test_archive as ota
test_data = ota.TestData()

'''
This test suite does not test command line usage. See
http://hg.python.org/cpython/file/default/Lib/test/test_cmd_line_script.py for
several examples of how to test handling of command line arguments.
'''

class SARWindTest(unittest.TestCase):
    def setUp(self):
        if test_data.noData:
            raise ValueError('No test data available')

    def test_sarwind_using_default(self):
        for i in range(len(test_data.radarsat2)):
            w = SARWind(test_data.radarsat2[i])
            self.assertIsInstance(w, SARWind)

    def test_sarwind_using_filenames(self):
        for i in range(len(test_data.asar)):
            w = SARWind(test_data.asar[i], wind_direction=test_data.ncep4asar[i])
            w.plot(filename='asar_agulhas_test_plot.png', show=False)
            self.assertIsInstance(w, SARWind)

    def test_sarwind_using_asar_filename_ncep_nansat(self):
        for i in range(len(test_data.asar)):
            mw = Nansat(test_data.ncep4asar[i])
            w = SARWind(test_data.asar[i], wind_direction=mw)
            self.assertIsInstance(w, SARWind)

    def test_sarwind_using_asar_nansat_ncep_filename(self):
        for i in range(len(test_data.asar)):
            asar = Nansat(test_data.asar[i])
            w = SARWind(asar, wind_direction=test_data.ncep4asar[i])
            self.assertIsInstance(w, SARWind)

    def test_sarwind_using_asar_nansat_ncep_nansat(self):
        for i in range(len(test_data.asar)):
            asar = Nansat(test_data.asar[i])
            mw = Nansat(test_data.ncep4asar[i])
            w = SARWind(asar, wind_direction=mw)
            self.assertIsInstance(w, SARWind)

    def test_plot_agulhas(self):
        '''
        Tests the SARWind.plot method. Remember to update the file
        basic_agulhas_ref.png if plot method is changed.

        If adding more test plots, please make sure the indices are correct
        '''
        w = SARWind(test_data.asar[0], wind_direction=test_data.ncep4asar[0])
        w.plot(filename=os.path.join(dirname,
            'plots/agulhas_test_plot.png'), show=False)
        self.assertTrue(filecmp.cmp(
            os.path.join(dirname,'plots/agulhas_test_plot.png'),
            os.path.join(dirname,'plots/agulhas_test_plot_ref.png')))

    def test_plot_barents(self):
        w = SARWind(test_data.radarsat2[0])
        w.plot(filename=os.path.join(dirname,
            'plots/rs2_barents_test_plot.png'), show=False)
        self.assertTrue(filecmp.cmp(
            os.path.join(dirname,'plots/rs2_barents_test_plot.png'),
            os.path.join(dirname,'plots/rs2_barents_test_plot_ref.png')))

    def tearDown(self):
        # Delete test plots
        if os.path.exists( os.path.join( dirname,
                'plots/agulhas_test_plot.png' ) ):
            os.unlink(os.path.join(dirname,'plots/agulhas_test_plot.png'))
        if os.path.exists( os.path.join( dirname,
                'plots/rs2_barents_test_plot.png' ) ):
            os.unlink(os.path.join(dirname,'plots/rs2_barents_test_plot.png'))


if __name__=='__main__':
    unittest.main()
