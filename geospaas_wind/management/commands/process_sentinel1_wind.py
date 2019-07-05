import pytz
from dateutil.parser import parse
from django.utils import timezone
from django.core.management.base import BaseCommand, CommandError
from django.core.exceptions import MultipleObjectsReturned

from geospaas.catalog.models import Dataset

from geospaas_wind.models import Dataset as WindDataset

class Command(BaseCommand):
    help = """
        Process high resolution wind from Sentinel-1 SAR availabe on OPeNDAP

        Args:
            --date <date>: Select Sentinel-1 datasets by date (yyyy-mm-dd)
            --url <str>: Select Sentinel-1 dataset by url
            --force-reprocessing: Force reprocessing of SAR wind

        Example: 
            Process all Sentinel-1 datasets on 2019-07-01 (the datasets must have been added to the
            Geo-SPaaS catalog first):

            ./manage.py process_sentinel1_wind --date 2019-07-01

        """

    def add_arguments(self, parser):
        parser.add_argument('--date', type=str, default='', 
                help='Date of coverage (yyyy-mm-dd)')
        parser.add_argument('--force-reprocessing', action='store_true', default=False,
                help='Force reprocessing')
        parser.add_argument('--url', type=str, default='', help='URL of OPeNDAP dataset')


    def handle(self, *args, **options):
        if options['url']:
            s1ds = Dataset.objects.filter(dataseturi__uri=options['url'])
        else:
            s1ds = Dataset.objects.filter(source__platform__short_name__contains='SENTINEL-1')
        if options['date']:
            tt = parse(options['date'])
            tt = timezone.datetime(tt.year, tt.month, tt.day, tzinfo=pytz.utc)
            s1ds = s1ds.filter(time_coverage_start__range=[tt, tt + timezone.timedelta(hours=24)])

        failed = 0
        reprocessed = 0
        ingested = 0
        numds = len(s1ds)
        print('Processing %d datasets' %numds)
        for i,ds in enumerate(s1ds):
            try:
                dsuri = ds.dataseturi_set.get(service='OPENDAP')
            except MultipleObjectsReturned:
                failed += 1
                continue
            wds, cr = WindDataset.objects.process(dsuri.uri, force=options['force_reprocessing'])
            if cr:
                self.stdout.write('Successfully processed (%d/%d): %s\n' % (i+1, numds, dsuri.uri))
                ingested += 1
            elif not cr and wds:
                self.stdout.write('Successfully reprocessed (%d/%d): %s\n' % (i+1, numds, dsuri.uri))
                reprocessed += 1
            elif not cr and not wds:
                failed += 1

        self.stdout.write('Successfully processed %d/%d datasets.\n' % (ingested, numds))
        if  options['force_reprocessing']:
            self.stdout.write('Successfully reprocessed %d/%d datasets.\n' % (reprocessed, numds))
        self.stdout.write('Failed at processing %d/%d datasets.\n' % (failed, numds))
