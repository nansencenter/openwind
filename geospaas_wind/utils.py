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
import urllib.request

from datetime import datetime, timedelta
from dateutil.parser import parse

from nansat.nansat import Nansat

from openwind.sar_wind import SARWind

from django.core.management import call_command
from django.contrib.gis.geos import WKTReader

from geospaas.utils.utils import nansat_filename
from geospaas.nansat_ingestor.models import Dataset as NansatDataset
from geospaas.catalog.models import Dataset


# Place to store downloads - TODO: generalise to follow env settings
downloads = os.path.join(os.path.expanduser('~'), 'downloads')
if not os.path.isdir(downloads):
    os.mkdir(downloads)

def wind_from_downloaded_sar_and_arome_netcdfs(sar, arome):
    if type(sar)==str:
        sar = Dataset.objects.get(dataseturi__uri=sar)
    if type(arome)==str:
        arome = Dataset.objects.get(dataseturi__uri=arome)
    # Get url of downloadable s1 netcdf
    sar_nc_uri = sar.dataseturi_set.get(service='HTTPServer')
    sar_filename = os.path.join(downloads, os.path.basename(sar_nc_uri.uri))
    # Download the netcdf for temporary storage
    urllib.request.urlretrieve(sar_nc_uri.uri, sar_filename)

    # Get url of downloadable arome netcdf
    uri = arome.dataseturi_set.get(service='HTTPServer')
    arome_fn = os.path.join(downloads, os.path.basename(uri.uri))
    # Download the netcdf for temporary storage
    urllib.request.urlretrieve(uri.uri, arome_fn)

    # Calculate wind
    w = SARWind(sar_filename, arome_fn)

    # Remove the wind field data from cache
    os.remove(arome_fn)
    os.remove(sar_filename)

    return w

def wind_from_sar_and_arome_forecast(sar_uri):
    """ Calculate wind field from provided SAR and forecast wind files
    """
    result_fn = os.path.join(downloads, 'wind_from_'+os.path.basename(sar_uri))
    sar_ds = NansatDataset.objects.get(dataseturi__uri=sar_uri)
    # Find collocated Arome forecast
    arome_ds = Dataset.objects.filter(
            summary__contains='AROME', 
            time_coverage_start__range=[
                sar_ds.time_coverage_start - timedelta(hours=3), 
                sar_ds.time_coverage_end + timedelta(hours=3)],
            geographic_location__geometry__intersects=sar_ds.geographic_location.geometry
        ).order_by('time_coverage_start')[0]

    try:
        # Calculate wind
        w = SARWind(sar_uri, arome_ds.dataseturi_set.get(service='OPENDAP'))
    except:
        """ Likely a problem with the server - fall back to downloading the data """
        w = wind_from_downloaded_sar_and_arome_netcdfs(sar_ds, arome_ds)

    # Export where the data is used... Not always useful to store netCDFs
    #w.export(result_fn)
    return w

def crawl_arome_arctic_archive(date, base_url='https://thredds.met.no/thredds/catalog/aromearcticarchive/'):
    """ Crawl met.no thredds server and ingest datasets for given date into the catalog
    """
    url = base_url + date.strftime('%Y/%m/%d') + '/catalog.html'
    call_command('ingest_thredds_crawl', url, filename='arome_arctic_vtk*')

def crawl_sentinel1_archive(date, base_url='http://nbstds.met.no/thredds/catalog/NBS/'):
    """ Crawl met.no thredds server and ingest datasets for given date into the catalog
    """
    urlS1A = base_url + 'S1A/' + date.strftime('%Y/%m/%d') + '/catalog.html'
    urlS2A = base_url + 'S1B/' + date.strftime('%Y/%m/%d') + '/catalog.html'
    call_command('ingest_thredds_crawl', urlS1A)
    call_command('ingest_thredds_crawl', urlS2A)

def crawl_data_archives(date):
    crawl_arome_arctic_archive(date)
    crawl_sentinel1_archive(date)
