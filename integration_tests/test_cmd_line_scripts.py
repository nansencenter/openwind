#-------------------------------------------------------------------------------
# Name:
# Purpose:      
#
# Author:       Morten Wergeland Hansen
# Modified:	Morten Wergeland Hansen
#
# Created:	07.04.2015
# Last modified:07.04.2015 13:51
# Copyright:    (c) NERSC
# License:      
#-------------------------------------------------------------------------------
import unittest

from openwind.sar_wind import create_parser, SARWind
import openwind.tests.openwind_test_archive as ota

'''
See http://dustinrcollins.com/testing-python-command-line-apps
'''

class CmdLineTests(unittest.TestCase):

    @classmethod
    def setUpClass(cls):
        parser = create_parser()
        cls.parser = parser
        cls.test_data = ota.OpenWindTestData()
        if cls.test_data.noData:
            raise ValueError('No test data available')

    def test_with_empty_args(self):
        '''
        User passes no args, should fail with SystemExit
        '''
        with self.assertRaises(SystemExit):
            self.parser.parse_args([])

    def test_constant_wind_dir(self):
        '''
        User passes constant wind direction
        '''
        if len(self.test_data.radarsat2)==0:
            raise IOError('No Radarsat-2 data - try adding some as ' \
                    'described in templates/openwind_local_archive.py' )
        args = self.parser.parse_args([self.test_data.radarsat2[0], '-f',
            'test.jpg', '-w', '90'])
        w = SARWind(args.sar_image, args.wind_direction, pixelsize=args.pixelsize)
        self.assertIsInstance(w, SARWind)


if __name__=='__main__':
    unittest.main()
