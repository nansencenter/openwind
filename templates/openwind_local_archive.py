#-------------------------------------------------------------------------------
# Name:         openwind_local_archive.py
# Purpose:      Add the paths to local test data in this file and place it on
#               your python path.
#
# Author:       Morten Wergeland Hansen
# Modified:	Morten Wergeland Hansen
#
# Created:	26.05.2014
# Last modified:16.06.2014 09:46
# Copyright:    (c) NERSC
# License:      
#-------------------------------------------------------------------------------
import os
path = os.path.dirname(__file__)

# Radarsat-2 file:
rs2 = '/Volumes/sat_downloads_radarsat2/RS2_20140326_054601_0076_SCWA_HHHV_SGF_316922_8221_9516957.zip'
try:
    with open(rs2) as f:
        pass
except IOError as exc:
    rs2 = None
    print "I/O Error (%s): %s" %(rs2, exc.strerror)
    #raise IOError("%s: %s" % (rs2, exc.strerror))

# Quad-pol Radarsat-2
rs2_quad = '/Volumes/sat_downloads_radarsat2/fine_quad_pol/RS2_20131214_084140_0004_FQ4_HHVVHVVH_SLC_299782_1914_9217833.zip'
try:
    with open(rs2_quad) as f:
        pass
except IOError as exc:
    rs2_quad = None
    print "I/O Error (%s): %s" %(rs2, exc.strerror)
