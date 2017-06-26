# Test Sentinel-1 SAR

import unittest

from openwind.sar_wind import SARWind
from nansat.nansat import Nansat, Domain

import openwind_integration_tests.openwind_test_archive as ota

class S1Test(unittest.TestCase):

    def setUp(self):
        #self.test_data = ota.OpenWindTestData()
        #self.test_data.get_sentinel1a_small()
        self.s1aEW = '/vagrant/shared/test_data/sentinel1_l1/S1A_EW_GRDM_1SDH_20170227T065537_20170227T065637_015466_019652_DD9C.SAFE'
        self.s1bIW = '/vagrant/shared/test_data/sentinel1_l1/S1B_IW_GRDM_1SDV_20170227T061040_20170227T061113_004482_007CD7_3205.SAFE'
        self.arome_arctic = '/vagrant/shared/test_data/generic/arome_arctic_pp_2_5km_20170227T00Z.nc'
        self.arome_metcoop = '/vagrant/shared/test_data/generic/arome_metcoop_default2_5km_20170227T00Z.nc'
        self.ecmwf_metno = '/vagrant/shared/test_data/generic/ec_atmo_0_1deg_20170227T000000Z_1h.nc'

    def test_s1a_without_explicit_wind_direction(self):
        w = SARWind(self.s1bIW)
        self.assertIsInstance(w, SARWind)

    def test_s1aEW_with_arome_arctic(self):
        w = SARWind(self.s1aEW, wind_direction=self.arome_arctic)
        self.assertIsInstance(w, SARWind)

    def test_s1bIW_with_arome_metcoop(self):
        w = SARWind(self.s1aEW, wind_direction=self.arome_metcoop)
        self.assertIsInstance(w, SARWind)

    def test_s1aEW_with_ecmwf(self):
        w = SARWind(self.s1aEW, wind_direction=self.ecmwf_metno)
        self.assertIsInstance(w, SARWind)

    def test_s1bIW_with_ecmwf(self):
        w = SARWind(self.s1aEW, wind_direction=self.ecmwf_metno)
        self.assertIsInstance(w, SARWind)
