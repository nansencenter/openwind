import os, warnings
import json
import numpy as np
import tempfile

import pythesint as pti

from nansat.nansat import Nansat
from nansat.domain import Domain
from nansat.tools import haversine

from django.db import models
from django.conf import settings
from django.utils import timezone
from django.core.exceptions import ObjectDoesNotExist

from geospaas.catalog.models import DatasetURI, Dataset
from geospaas.nansat_ingestor.managers import DatasetManager

from geospaas_wind.utils import wind_from_sar_and_arome_forecast
from geospaas_wind.exceptions import TooHighResolutionError, PolarizationError

class WindManager(DatasetManager):

    def process(self, uri, force=False, *args, **kwargs):
        fn = 'WIND_'+os.path.basename(uri)
        if DatasetURI.objects.filter(uri__contains=fn) and not force:
            wds = Dataset.objects.filter(dataseturi__uri__contains=fn)[0]
            return wds, False

        thredds_fn = os.path.join(settings.PRODUCTS_ROOT, fn)
        wind_uri = 'file://localhost' + thredds_fn

        try:
            w = wind_from_sar_and_arome_forecast(uri)
        except (TooHighResolutionError, PolarizationError, ObjectDoesNotExist) as e:
            if type(e)==Dataset.DoesNotExist:
                warnings.warn(uri+': '+e.args[0])
            else:
                # ObjectDoesNotExist could happen if there is no overlap between SAR and model
                warnings.warn(e.file + ': ' + e.msg)
            return None, False

        metadata = w.get_metadata()

        # Direct reprojection fails - gdal can't read the bands if we do w.reproject...
        # Workaround: Export wind to temporary file
        #tmp_filename = os.path.join(settings.PRODUCTS_ROOT,'TMP_WIND_'+os.path.basename(uri))
        fd, tmp_filename = tempfile.mkstemp(suffix='.nc')
        os.close(fd) # Just in case - see https://www.logilab.org/blogentry/17873
        w.export(tmp_filename)

        # Read temporary file
        ww = Nansat(tmp_filename)

        # Reproject
        lon, lat = ww.get_geolocation_grids()
        #lon, lat = w.get_geolocation_grids()
        srs = '+proj=stere +datum=WGS84 +ellps=WGS84 +lat_0=%.2f +lon_0=%.2f +no_defs'%(np.mean(lat),np.mean(lon))
        xmin, xmax, ymin, ymax = -haversine(np.mean(lon),np.mean(lat),np.min(lon),np.mean(lat)), \
                    haversine(np.mean(lon),np.mean(lat),np.max(lon),np.mean(lat)), \
                    -haversine(np.mean(lon),np.mean(lat),np.mean(lon),np.min(lat)), \
                    haversine(np.mean(lon),np.mean(lat),np.mean(lon),np.max(lat))
        ext = '-te %d %d %d %d -tr 500 500' %(xmin-10000, ymin-10000, xmax+10000, ymax+10000)
        d = Domain(srs, ext)
        ww.reproject(d, tps=True)
        #w.reproject(d, tps=True)

        # Set global metadata
        metadata['data_center'] = json.dumps(pti.get_gcmd_provider('nersc'))
        metadata['naming_authority'] = 'no.nersc.sios_infranor'
        metadata['project'] = 'SIOS InfraNor'
        metadata['entry_title'] = 'Wind field from '+os.path.basename(uri)
        metadata.pop('file_creation_date')
        metadata['history'] = metadata['history'] + ' ' + timezone.now().isoformat() + \
                '. Calculated wind field from NRCS and Arome Arctic forecast wind directions.'
        metadata.pop('institution')
        metadata['keywords'] += ', ['
        for key, value in pti.get_gcmd_science_keyword('U/V WIND COMPONENTS').items():
            if value:
                metadata['keywords'] += value + ', '
        metadata['keywords'] += ']'
        metadata.pop('LINE_SPACING')
        metadata.pop('PIXEL_SPACING')
        metadata['summary'] = 'Near surface (10m) wind from Arome Arctic forecast wind and ' + metadata['summary']
        metadata['title'] = 'Near surface wind from '+metadata['title']

        # Export
        ww.export2thredds(thredds_fn, mask_name='swathmask', metadata=metadata, no_mask_value=1)
        wds, cr = super(WindManager, self).get_or_create(wind_uri)
        
        os.unlink(tmp_filename)

        return wds, cr

