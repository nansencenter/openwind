import json
import numpy as np
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

        # Remove GCPs 
        #self.dataset.SetGCPs([],'')
        # Set projection to wkt
        self.dataset.SetProjection(NSR().wkt)


        band_lat = self.dataset.GetRasterBand(1)
        # Check that it is actually longitudes
        if not band_lat.GetMetadata()['standard_name'] == 'latitude':
            raise ValueError('Cannot find latitude band')
        lat = band_lat.ReadAsArray()

        band_lon = self.dataset.GetRasterBand(2)
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
        src = {'SourceFilename': self.band_vrts['new_lon_VRT'].filename,
               'SourceBand': 1}
        dst = {'wkv': 'longitude',
               'name': 'lon_corrected'}

        # Find latitude band number
        fn = self.sub_filenames(args[1])
        lat_band_num = [ii for ii, ll in enumerate(fn) if ':lat' in ll][0] + 1

        self.dataset.SetProjection(NSR().wkt)
        #self.dataset.SetGeoTransform((0.0, 1.0, 0.0, 0.0, 0.0, 1.0))

        self.dataset.SetGCPs(VRT._lonlat2gcps(lon, lat, n_gcps=400), NSR().wkt)

        # Add geolocation from correct longitudes and latitudes
        self._add_geolocation(Geolocation(self.band_vrts['new_lon_VRT'], self, 1, lat_band_num)) # band numbers are hardcoded...
        

        # TODO: add GCMD/DIF metadata

    def _create_empty(self, gdal_dataset, gdal_metadata):
        fn = self.sub_filenames(gdal_dataset)
        if not fn:
            raise WrongMapperError
        #latfn = [ll for ll in fn if ':lat' in ll][0]
        #lonfn = [ll for ll in fn if ':lon' in ll][0]
        # This may destroy things..
        #lon = np.mod(gdal.Open(lonfn).ReadAsArray() + 180, 360) - 180
        #lon = gdal.Open(lonfn).ReadAsArray()
        #lat = gdal.Open(latfn).ReadAsArray()
        #self._init_from_lonlat(lon, lat)

        lat = gdal.Open(fn[0])

        super(NetcdfCF, self).__init__(lat.RasterXSize, lat.RasterYSize, metadata=gdal_metadata)

        
