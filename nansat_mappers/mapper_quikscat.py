import json
import numpy as np
from datetime import datetime
import pythesint as pti

import gdal

from nansat.vrt import VRT
from nansat.geolocation import Geolocation
from nansat.nsr import NSR
from nansat.domain import Domain
from nansat.mappers.mapper_netcdf_cf import Mapper as NetcdfCF
from nansat.exceptions import WrongMapperError

class Mapper(NetcdfCF):

    def __init__(self, *args, **kwargs):
        quartile = kwargs.pop('quartile', 0)

        metadata = args[2]
        if not metadata.has_key('NC_GLOBAL#source') \
                or not metadata['NC_GLOBAL#source'].lower() == 'quikscat':
            raise WrongMapperError

        super(Mapper, self).__init__(*args, **kwargs)
        
        intervals = [0,1,2,3]
        if not quartile in intervals:
            raise ValueError('quartile must be one of [0,1,2,3]')

        y_size = self.dataset.RasterYSize/4
        y_offset = [y_size*qq for qq in intervals][quartile]

        # Crop
        self.set_offset_size('y', y_offset, y_size)

        # Create band of times
        # TODO: resolve nansat issue #263 (https://github.com/nansencenter/nansat/issues/263)
        tt = self.times()[y_offset : y_offset + y_size]
        self.dataset.SetMetadataItem('time_coverage_start', tt[0].astype(datetime).isoformat())
        self.dataset.SetMetadataItem('time_coverage_end', tt[-1].astype(datetime).isoformat())
        time_stamps = (tt - tt[0]) / np.timedelta64(1, 's')
        self.band_vrts['time'] = VRT.from_array(
                np.tile(time_stamps, (self.dataset.RasterXSize, 1)).transpose()
            )
        self.create_band(
                src = {
                    'SourceFilename': self.band_vrts['time'].filename,
                    'SourceBand': 1,
                },
                dst = {
                    'name': 'timestamp', 
                    'time_coverage_start': tt[0].astype(datetime).isoformat(), 
                    'units': 'seconds since time_coverage_start',
                }
            )

        # Set projection to wkt
        self.dataset.SetProjection(NSR().wkt)

        band_lat = self.dataset.GetRasterBand(self._latitude_band_number())
        # Check that it is actually latitudes
        if not band_lat.GetMetadata()['standard_name'] == 'latitude':
            raise ValueError('Cannot find latitude band')
        lat = band_lat.ReadAsArray()

        band_lon = self.dataset.GetRasterBand(self._longitude_band_number())
        # Check that it is actually longitudes
        if not band_lon.GetMetadata()['standard_name'] == 'longitude':
            raise ValueError('Cannot find longitude band')
        lon = band_lon.ReadAsArray()

        # Apply correction of longitudes (they are defined on 0:360 degrees but also contain
        # negative values)
        # TODO: consider making this core to nansat - different ways of defining longitudes
        # (-180:180 og 0:360 degrees) often cause problems...
        lon = np.mod(lon+180., 360.) - 180.

        self.band_vrts['new_lon_VRT'] = VRT.from_array(lon)

        self.dataset.SetProjection(NSR().wkt)

        self.dataset.SetGCPs(VRT._lonlat2gcps(lon, lat, n_gcps=400), NSR().wkt)

        # Add geolocation from correct longitudes and latitudes
        self._add_geolocation(Geolocation(self.band_vrts['new_lon_VRT'], self, 1,
            self._latitude_band_number())) # band numbers are hardcoded...

        # TODO: add GCMD/DIF metadata

    def _latitude_band_number(self):
        return [ii for ii, ll in enumerate(self.sub_filenames()) if ':lat' in ll][0] + 1

    def _longitude_band_number(self):
        return [ii for ii, ll in enumerate(self.sub_filenames()) if ':lon' in ll][0] + 1

    def _create_empty(self, gdal_dataset, gdal_metadata):
        lat = gdal.Open(self.sub_filenames()[self._latitude_band_number()])
        super(NetcdfCF, self).__init__(lat.RasterXSize, lat.RasterYSize, metadata=gdal_metadata)

        
