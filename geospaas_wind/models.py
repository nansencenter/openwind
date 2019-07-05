from django.db import models
from geospaas.catalog.models import Dataset as CatalogDataset

from some_new_package.managers import WindManager

class Dataset(CatalogDataset):
    class Meta:
        proxy = True

    objects = WindManager()
