from django.db import models
from geospaas.catalog.models import Dataset as CatalogDataset

from geospaas_wind.managers import WindManager

class Dataset(CatalogDataset):
    class Meta:
        proxy = True

    objects = WindManager()
