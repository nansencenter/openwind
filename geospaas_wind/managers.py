import os, warnings

from django.conf import settings
from django.db import models

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
        w.export2thredds(thredds_fn)
        wds, cr = super(WindManager, self).get_or_create(wind_uri)
        
        return wds, cr

