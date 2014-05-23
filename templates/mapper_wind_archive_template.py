# Name:         mapper_wind_archive.py
# Purpose:      Nansat mapping for local wind files
# Author:       Knut-Frode Dagestad
# Licence:      This file is part of OpenWind. You can redistribute it or
#               modify under the terms of GNU General Public License, v.3
#               http://www.gnu.org/licenses/gpl-3.0.html
#
# This is a template mapper. To adapt to your local computer:
# - copy this file to any folder within your PYTHONPATH
# - rename to any filename starting with 'mapper_" and ending with ".py"
# - Modify code below to find windFileName at your local file system
#
#
# Mapper searches local file archive for the Nansat-readable wind
# file closest to the requested time, specified as keyword as show below
#
# Usage:
#    w = Nansat('wind_archive:YYYYMMDDHHMM')
# Nansat example:
#    w = Nansat('wind_archive:201405011000')
# OpenWind example:
#    sw = SARWind(sar_file, winddir='wind_archive:201405011000')
#
# From commandline:
# $ nansatinfo wind_archive:201405011000
# $ sar_wind.py -s <sar_filename> -w wind_archive:201405011000 -f fig.png

from datetime import datetime, timedelta

from nansat.vrt import VRT
from nansat import Nansat

localFolder = '/disk2/data/modelwind/ncep/gfs/'
keywordBase = 'wind_archive'


class Mapper(VRT):

    def __init__(self, fileName, gdalDataset, gdalMetadata, **kwargs):
        ''' Create VRT '''

        ##############
        # Get time
        ##############
        if fileName[0:len(keywordBase)] != keywordBase:
            raise AttributeError("Wrong mapper")

        timestr = fileName[len(keywordBase)+1::]
        time = datetime.strptime(timestr, '%Y%m%d%H%M')

        ######################################################
        # Find windFileName corresponding to a Nansat-readable
        # file in your local (or remote) file archive
        ######################################################
        windFileName = localFolder + <.......>

        ######################################################
        # Open file with any other Nansat mapper
        ######################################################
        w = Nansat(windFileName)
        VRT.__init__(self, vrtDataset=w.vrt.dataset)

        return
