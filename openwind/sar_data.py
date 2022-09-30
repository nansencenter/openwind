from openwind.utils import measure_time, check_inputs
import zipfile
import asf_search as asf
from s1denoise import Sentinel1Image
from nansat import Nansat
from pathlib import Path
import numpy as np
from typing import Union, Optional


def _retrieve_asf_creds():
    with open('/home/artmoi/.asfapirc') as f:
        lines = [line.strip().split(': ') for line in f.readlines()]
        creds = {line[0]:line[1] for line in lines}
    return creds


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
        denoise_alg: Optional[str] = None, 
        dst_px_size: Optional[float] = None
    ) -> Nansat:
    """
    Read, denoise, and/or resize SAR data (tested for S1, RS2, and ASAR)
    NOTE: Additional denoising scheme is implemented only for Sentinel-1 data

    :param sar_product_uri: path/to/sar/product
    :param denoise_alg: name of denoising scheme. 
        NOTE: Only rquired for Sentinel-1 data
    :param dst_px_size: target pixel size in meters
    :returns sar_data: sar dataset opened in nansat
    """
    print(f'>> Processing {sar_product_uri}')
    # Sentinel-1 requires additional denoising and hence must be threated separately
    if sar_product_uri.name.startswith('S1') and denoise_alg is not None:
        # Read the Sentinel-1 GRD file and get all aux calibration data
        sar_data = Sentinel1Image(str(sar_product_uri), mapperName='sentinel1_l1')
        # Remove thermal noise from the sigma0 and add it as a separate band to the dataset
        print(f'>> Remove noise from sigma0')
        sar_data.add_band(sar_data.remove_thermal_noise('VV', algorithm=denoise_alg),
                          parameters={'name':'sigma0_vv_denoised', 'algorithm': denoise_alg})
    # In case of other SAR data supported by nansat (e.g., RS2, ASAR) no additional 
    # calibrations applied. 
    else:
        sar_data = Nansat(str(sar_product_uri))
    # Resize (increase/decrease px size) image
    if dst_px_size is not None:
        print(f'>> Resizing image to {dst_px_size} m')
        # Calculate resize factor based on source px size and target px size
        resize_factor = np.mean(sar_data.get_pixelsize_meters()) / dst_px_size
        # Resize image to target resolution using bilinear resampling
        sar_data.resize(resize_factor, resample_alg=1)
    
    return sar_data
