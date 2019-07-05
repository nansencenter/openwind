from django.db import models

class WindManager(models.Manager):

    def get_or_create(self, sar_uri, wind_uri, *args, **kwargs):

