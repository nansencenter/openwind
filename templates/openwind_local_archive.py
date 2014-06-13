#-------------------------------------------------------------------------------
# Name:         openwind_local_archive.py
# Purpose:      Add the paths to local test data in this file and place it on
#               your python path.
#
# Author:       Morten Wergeland Hansen
# Modified:	Morten Wergeland Hansen
#
# Created:	26.05.2014
# Last modified:12.06.2014 15:59
# Copyright:    (c) NERSC
# License:      
#-------------------------------------------------------------------------------
import os
path = os.path.dirname(__file__)

# A Radarsat-2 file:
rs2 = '/Volumes/sat_downloads_radarsat2/RS2_OK53265_PK506736_DK449935_SCWA_20140610_134945_HH_SGF.ZIP'
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
