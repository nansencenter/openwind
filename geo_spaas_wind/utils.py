""" Utility functions for finding SAR and auxiliary wind information data to process SAR wind

All data management, i.e., finding and collocating data, should be done here.

The following options for auxiliary wind direction information are covered:

    1) constant wind direction - use sar_wind_from_constant
    2) the name of a Nansat compatible file containing wind direction information
        - Arome arctic
        - Arome metcoop
        - NCEP
        - HIRLAM
    3) using nansat.mappers.mapper_ncep_wind_online

The SAR dataset is provided as either a filename or a Nansat object.

The wind field is calculated directly from CMOD, and optionally also by using Bayesian statistics.

"""
import os
from datetime import datetime, timedelta
from dateutil.parser import parse

from nansat.nansat import Nansat

from openwind.sar_wind import SARWind

from django.core.management import call_command
from django.contrib.gis.geos import WKTReader

from geospaas.catalog.models import Dataset


# Place to store downloads - TODO: generalise to follow env settings
downloads = os.path.join(os.path.expanduser('~'), 'downloads')
if not os.path.isdir(downloads):
    os.mkdir(downloads)

def wind_from_sar_and_arome_forecast(sar_filename):
    """ Calculate wind field from provided SAR and forecast wind files
    """
    n = Nansat(sar_filename)
    geometry = WKTReader().read(n.get_border_wkt(n_points=1000))
    arome_ds = Dataset.objects.filter(
            summary__contains='AROME', 
            time_coverage_start__range=[
                parse(n.get_metadata('time_coverage_start'))-timedelta(hours=3), 
                parse(n.get_metadata('time_coverage_end'))+timedelta(hours=3)],
            geographic_location__geometry__intersects=geometry
        ).order_by('time_coverage_start')[0]
    import ipdb
    ipdb.set_trace()
    return SARWind(sar_filename, arome_ds.dataseturi_set.get().uri)

def sar_wind_from_ncep_online(sar_filename):
    """ Calculate wind field from provided SAR file and online ncep wind
    """
    raise NotImplementedError

def sar_wind_from_constant(sar_filename, wind_direction):
    """ Calculate wind field from provided SAR file and a constant wind direction relative to the
    SAR look direction
    """
    raise NotImplementedError

def crawl_arome_arctic_archive(date, base_url='https://thredds.met.no/thredds/catalog/aromearcticarchive/'):
    """ Crawl met.no thredds server and ingest datasets for given date into the catalog
    """
    url = base_url + date.strftime('%Y/%m/%d') + '/catalog.html'
    call_command('ingest_thredds_crawl', url, filename='arome_arctic_vtk*')
