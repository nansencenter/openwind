# Test Sentinel-1 SAR

import unittest

from openwind.sar_wind import SARWind
from nansat.nansat import Nansat, Domain

import openwind_integration_tests.openwind_test_archive as ota

class S1ATest(unittest.TestCase):

    def setUp(self):
        self.test_data = ota.OpenWindTestData()
        self.test_data.get_sentinel1a_small()

    def test_s1a_small(self):
        w = SARWind(self.test_data.sentinel1a['small'])
        self.assertIsInstance(w, SARWind)
