import os, warnings
import numpy as np

from django.conf import settings
from django.db import models

from nansat.domain import Domain
from nansat.tools import haversine

from geospaas.catalog.models import DatasetURI, Dataset

from geospaas_wind.utils import wind_from_sar_and_arome_forecast
from geospaas_wind.exceptions import TooHighResolutionError, PolarizationError

class WindManager(models.Manager):

    def process(self, uri, force=False, *args, **kwargs):
        fn = 'WIND_'+os.path.basename(uri)
        if DatasetURI.objects.filter(uri__contains=fn) and not force:
            wds = Dataset.objects.filter(dataseturi__uri__contains=fn)[0]
            return wds, False

        thredds_fn = os.path.join(settings.PRODUCTS_ROOT, fn)
        wind_uri = 'file://localhost' + thredds_fn

        try:
            w = wind_from_sar_and_arome_forecast(uri)
        except (TooHighResolutionError, PolarizationError) as e:
            warnings.warn(e.file + ': ' + e.msg)
            return None, False
        # Reproject
        lon, lat = w.get_geolocation_grids()
        srs = '+proj=stere +datum=WGS84 +ellps=WGS84 +lat_0=%.2f +lon_0=%.2f +no_defs'%(np.mean(lat),np.mean(lon))
        xmin, xmax, ymin, ymax = -haversine(np.mean(lon),np.mean(lat),np.min(lon),np.mean(lat)), \
                    haversine(np.mean(lon),np.mean(lat),np.max(lon),np.mean(lat)), \
                    -haversine(np.mean(lon),np.mean(lat),np.mean(lon),np.min(lat)), \
                    haversine(np.mean(lon),np.mean(lat),np.mean(lon),np.max(lat))
        ext = '-te %.2f %2.f %.2f %.2f -tr 500 500' %(xmin,ymin,xmax,ymax)
        d = Domain(srs, ext)
        ## Copy to new object with only a few bands? -- This is done in export2thredds..
        #    super(SARWind, self).from_domain(sar_image, *args, **kwargs)
        #    self.vrt = sar_image.vrt
        #    self.mapper = sar_image.mapper
        #    self.logger = sar_image.logger
        w.reproject(d)
        # Export
        w.export2thredds(thredds_fn, bands={'U':{},'V':{}})
        wds, cr = super(WindManager, self).get_or_create(wind_uri)
        
        return wds, cr

