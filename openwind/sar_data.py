import logging
import zipfile
from pathlib import Path
from s1denoise import Sentinel1Image
from s1denoise.tools import run_correction
from typing import Union, Optional

import asf_search as asf
import dateutil.parser
import numpy as np
import pythesint as pti
import xarray as xr
from nansat import Nansat
from numpy.typing import NDArray

from openwind.utils import measure_time, check_inputs, Extent


logger = logging.getLogger(__name__)


def _retrieve_asf_creds(src):
    with open(src) as f:
        lines = [line.strip().split(': ') for line in f.readlines()]
        creds = {line[0]:line[1] for line in lines}
    return creds


def db2linear(s0: NDArray) -> NDArray:
    """
    Convert sigma0 from dB to linear units

    :param s0:   sigma0 in dB
    :returns s0: sigma0 in linear units
    """
    s0 = 10**(s0/10)
    return s0


@measure_time
def _download_data(query_set, dst, asf_creds):
    # If asf creds are not provided then retrieve from the system file
    asf_creds = _retrieve_asf_creds() if asf_creds is None else asf_creds
    # Create a session with provided login and password
    session = asf.ASFSession().auth_with_creds(asf_creds['login'], asf_creds['password'])
    # Download data
    query_set.download(path=dst, session=session)


@measure_time
def _unzip_data(uri, rm_zip=False):
    print(f'>> Unpacking {uri}', end=' ')
    with zipfile.ZipFile(uri, 'r') as s1_zipfile:
        s1_zipfile.extractall(path=uri.parent)
    # Remove original zip file if required
    uri.unlink() if rm_zip else None
    # Return uri to SAFE files
    return uri.with_suffix('.SAFE')


@check_inputs
def fetch_asf_s1_data(granule_id=None, query=None, download=False, dst=None, asf_creds=None, unpack=False, force=False):
    """
    Find and download (optinally) Sentinel-1 products from the ASF (https://asf.alaska.edu/) database
    using asf_search API using query set or granule id (see https://docs.asf.alaska.edu/asf_search/searching/)

    :param granule_id: str
    :param query (dict): dictionary with search parameters for the asf API
    :param download (bool): if True downloads all found products to dst path
    :param dst (str): 'path/to/downloads/dir'
    :param asf_creds (dict): dictionary with login and password to access the ASF follow
    :param unpack (bool): each product is a zip archive that will be unpacked to SAFE directory, all
        original zip files will be removed
    :returns query_set: list of products from ASF request or list or uris if download is True
    """

    if granule_id is not None:
        query_set = asf.granule_search(granule_id)
    # If granule id is not provided then search data using query parameters
    elif query is not None:
        query_set = asf.geo_search(**query)

    # If no data found then raise error
    if len(query_set) == 0 or download == False:
        uris = list()

    print(f">> Found {len(query_set)} products in ASF database")

    if download:
        # ensure that the dst is a path object
        dst = Path(dst)
        print(f'>> Downloading {len(query_set)} files:', end=' ')
        _download_data(query_set, dst, asf_creds)
        # Generate a list of downloaded uris
        uris = [dst / s1_product.properties['fileName'] for s1_product in query_set]

    if unpack:
        # Unzip each file and return list of SAFE uris
        uris = [_unzip_data(uri) for uri in uris]

    return query_set, uris


def preprocess_sar_data(
        sar_product_uri: Path,
        denoise_alg: str = 'NERSC',
        dst_px_size: Optional[float] = None
    ) -> Nansat:
    """
    Read, denoise, and/or resize SAR data (tested for S1, RS2, and ASAR)
    NOTE: Additional denoising scheme is implemented only for Sentinel-1 data

    :param sar_product_uri: path/to/sar/product
    :param denoise_alg: name of denoising scheme.
        NOTE: Only rquired for Sentinel-1 data
        NOTE: Currently only uses default NERSC algorithm
    :param dst_px_size: target pixel size in meters
    :returns sar_data: sar dataset opened in nansat
    """
    print(f'>> Processing {sar_product_uri}')
    # Sentinel-1 requires additional denoising and hence must be threated separately
    if sar_product_uri.name.startswith('S1') and denoise_alg is not None:
        # Read the Sentinel-1 GRD file and get all aux calibration data
        sar_data = Nansat(str(sar_product_uri), mapperName='sentinel1_l1')
        # Remove thermal noise from the sigma0 and add it as a separate band to the dataset
        # TODO: Find the way to specify correction algorithm: currently used default NERSC
        print(f'>> Remove noise from sigma0')

        denoised_s0 = run_correction(str(sar_product_uri))
        # TODO: dynamically specify polarization
        for pol in ['VV']: #denoised_s0:
            sar_data.add_band(denoised_s0[pol], parameters={
                'name': f'sigma0_{pol.lower()}_denoised',
                'algorithm': denoise_alg, 'units': 'dB'
            })
    # In case of other SAR data supported by nansat (e.g., RS2, ASAR) no additional
    # calibrations applied.
    else:
        sar_data = Nansat(str(sar_product_uri))
    # Resize (increase/decrease px size) image
    if dst_px_size is not None:
        print(f'>> Resizing image to {dst_px_size} m')
        # Calculate resize factor based on source px size and target px size
        resize_factor = np.mean(sar_data.get_pixelsize_meters()) / dst_px_size
        # Resize image to target resolution using neares neibour
        sar_data.resize(resize_factor, resample_alg=0)

    return sar_data


def _prep_band_metadata(src_name: str) -> dict:
    # Get CF metadata for the band
    meta_dict = pti.get_cf_standard_name(src_name)
    # Add FillValue
    meta_dict['_FillValue'] = -999.
    return meta_dict


def export2netcdf(
        sar_ds: Nansat,
        dst_path: Union[str, Path] = Path('.')
    ) -> xr.Dataset:
    """
    Export Nansat dataset to xarray an  d write NetCDF (Optional)

    :param sar_ds:   Nansat dataset
    :param dst_path: path/to/dst/file/dir
    """
    print(f'>> Fetch geolocation grids and watermask')
    # Extract geolocation grids for SAR frame
    lon_grd, lats_grd = sar_ds.get_geolocation_grids()
    # Extract MOD44W land/watermask for SAR frame
    watermask = sar_ds.watermask()
    # Create output xarray dataset
    print(f'>> Prepare dataset')
    out_ds = xr.Dataset(
        # Add metadata from original SAR acquisition
        attrs=sar_ds.get_metadata(),
        # Add coordinate bands
        coords=dict(
            lat=(('row', 'col'), lats_grd, _prep_band_metadata('latitude')),
            lon=(('row', 'col'), lon_grd, _prep_band_metadata('longitude'))),
        data_vars=dict(
            sigma0_s1dn=(('row', 'col'), db2linear(sar_ds['sigma0_vv_denoised']), {'polarization': 'VV', '_FillValue': -999.}),
            sigma0=(('row', 'col'), sar_ds['sigma0_VV'], {'polarization': 'VV', '_FillValue': -999.}),
            angle_of_incidence=(('row', 'col'), sar_ds['incidence_angle'], _prep_band_metadata('angle_of_incidence')),
            look_direction=(('row', 'col'), sar_ds['look_direction'], {'_FillValue': -999.}),
            watermask=(('row', 'col'), watermask[1], {'source': 'MOD44W', '_FillValue': -999.})))
    # Generate dst path for writing new dataset
    out_fname = Path(sar_ds.name).with_suffix('.nc')
    # Generate full path
    dst_path = Path(dst_path) / out_fname
    if dst_path.exists(): dst_path.unlink()
    # Write NetCDF dataset to the disk
    print(f'>> Write to: {dst_path}')
    out_ds.to_netcdf(dst_path)
    out_ds.close()
    return dst_path


class Sentinel1Source():
    """"""
    def __init__(self,
                 asf_product: asf.ASFProduct = None,
                 data_path:Union[str, Path] = None):
        self._asf_product = asf_product
        self.data_path = Path(data_path) if data_path is not None else None
        self.denoised_path = None

        self._bounding_box = None
        self._start_time = None
        self._end_time = None

    @classmethod
    def from_asf(cls, extent, time_start, time_end, query=None):
        """"""
        polygon = ("POLYGON(("
            f"{extent.west} {extent.south},"
            f"{extent.east} {extent.south},"
            f"{extent.east} {extent.north},"
            f"{extent.west} {extent.north},"
            f"{extent.west} {extent.south}))")

        if query is None:
            query = {
                'platform': asf.PLATFORM.SENTINEL1,
                'start': time_start,
                'end': time_end,
                'beamMode': asf.BEAMMODE.IW,
                'processingLevel': asf.PRODUCT_TYPE.GRD_HD,
                'intersectsWith': polygon,
            }
        query_set = asf.search(**query)
        return [cls(asf_product) for asf_product in query_set]

    @property
    def identifier(self):
        """"""
        if self._asf_product:
            return self._asf_product.properties['sceneName']
        elif self.data_path:
            return self.data_path.stem
        elif self.denoised_path:
            return self.denoised_path.stem

    @property
    def bounding_box(self):
        """"""
        if self._bounding_box is None:
            min_lon = 180.
            max_lon = -180.
            min_lat = 90.
            max_lat = -90.
            for point in self._asf_product.geometry['coordinates'][0]:
                if point[0] < min_lon:
                    min_lon = point[0]
                if point[0] > max_lon:
                    max_lon = point[0]
                if point[1] < min_lat:
                    min_lat = point[1]
                if point[1] > max_lat:
                    max_lat = point[1]
            self._bounding_box = Extent(min_lon, max_lon, min_lat, max_lat)
        return self._bounding_box

    @property
    def start_time(self):
        """"""
        if self._start_time is None:
            self._start_time = dateutil.parser.parse(self._asf_product.properties['startTime'])
        return self._start_time

    @property
    def end_time(self):
        """"""
        if self._end_time is None:
            self._end_time = dateutil.parser.parse(self._asf_product.properties['stopTime'])
        return self._end_time

    @property
    def properties(self):
        """"""
        return self._asf_product.properties

    def download(self, out_dir, unzip=False):
        """"""
        target = Path(out_dir, self._asf_product.properties['fileName'])
        logger.info("Downloading to %s", target)
        if not target.exists():
            self._asf_product.download(str(out_dir))
        else:
            logger.info("Did not download, destination already exists: %s", target)
        self.data_path = target
        if unzip:
            self.unzip(out_dir)
        return target

    def unzip(self, out_dir):
        """"""
        safe_name = f"{self.identifier}.SAFE/"
        safe_path = Path(out_dir, safe_name)
        if not safe_path.exists() and zipfile.is_zipfile(self.data_path):
            logger.info("Unzipping %s", self.data_path)
            with zipfile.ZipFile(self.data_path) as zip_file:
                if safe_name in zip_file.namelist():
                    zip_file.extractall(out_dir)
                    self.data_path = safe_path
                else:
                    raise RuntimeError(f"Not a Sentinel-1 SAFE archive: {self.data_path}")
        else:
            logger.info("Already unzipped: %s", self.data_path)
        return self.data_path

    def _run_s1_correction(self,
            angular_scale_copol=-0.2, angular_scale_crpol=-0.025, angular_offset=34.5,
            polarizations=('HH', 'VV', 'HV', 'VH'), algorithm='NERSC', dtype=np.float32,
            **kwargs):
        """Modified s1denoise.tools.run_correction to be able to specify polarizations
        """
        s1_image = Sentinel1Image(str(self.data_path))
        scale = {
            'HH': angular_scale_copol,
            'HV': angular_scale_crpol,
            'VH': angular_scale_crpol,
            'VV': angular_scale_copol,
        }
        denoised = {}
        for pol in polarizations:
            if pol in s1_image.pols:
                inc = s1_image.get_geolocation_full_size(pol, 'incidenceAngle')
                denoised[pol] = s1_image.remove_texture_noise(pol, algorithm=algorithm)
                denoised[pol] = (
                    10 * np.log10(denoised[pol]) - scale[pol] * (inc - angular_offset)
                ).astype(np.float32)
        return denoised

    def denoise(self, out_dir, algorithm='NERSC', polarizations=('VV',)):
        """"""
        output_path = out_dir / f'denoised_{self.identifier}.nc'
        logger.info("Denoising %s (%s) using %s algorithm",
                    self.identifier, ','.join(polarizations), algorithm)

        if output_path.exists():
            logger.info("Denoised file already exists: %s", output_path)
        else:
            denoised = self._run_s1_correction(polarizations=polarizations, algorithm=algorithm)

            logger.info(f'Creating denoised dataset at {output_path}...')
            s1_orig = Nansat(str(self.data_path))
            for pol in denoised:
                pol_low = pol.lower()
                s1_orig.add_band(denoised[pol], parameters={
                    'name': f'sigma0_{pol_low}_denoised',
                    'algorithm': algorithm,
                    'units': 'dB'
                })

            s1_orig.resize(pixelsize=500, resample_alg=0)
            lon_grd, lats_grd = s1_orig.get_geolocation_grids()
            watermask = s1_orig.watermask()

            denoised_vars = {}
            for pol in denoised:
                pol_low = pol.lower()
                band = f"sigma0_{pol_low}_denoised"
                denoised_vars[band] = (
                    ('row', 'col'),
                    db2linear(s1_orig[band]),
                    {'polarization': pol, '_FillValue': -999.})

            s1_dataset = xr.Dataset(
                data_vars={
                    "sigma0": (('row', 'col'), s1_orig['sigma0_VV'], {'polarization': 'VV', '_FillValue': -999.}),
                    "watermask": (('row', 'col'), watermask[1], {'source': 'MOD44W', '_FillValue': -999.}),
                    "incidence": (('row', 'col'), s1_orig['incidence_angle']),
                    "look_direction": (('row', 'col'), s1_orig['look_direction'], {'_FillValue': -999.}),
                    **denoised_vars,
                },
                coords={
                    "lon": (('row', 'col'), lon_grd),
                    "lat": (('row', 'col'), lats_grd)
                }
            )
            s1_dataset.to_netcdf(output_path)

        self.denoised_path = output_path
        return output_path
