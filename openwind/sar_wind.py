from numpy.typing import NDArray
from typing import Union, Optional
from nansat import Nansat
from openwind.sar_data import preprocess_s1_data
from openwind.wind_data import preprocess_era5_data, wind2sar_direction
from openwind.gmf.cmod5n import cmod5n_inverse
from pathlib import Path
import numpy as np


def derive_sar_wind(
        wind_source: Union[str, Path, Nansat],
        sar_source: Union[str, Path, Nansat],
        pixel_size_m: Union[int, float] = None,
        denoise_alg: Optional[str] = None,
        export: bool = False
    ) -> NDArray:
    """
    Derive wind product from collocated SAR observations and model field using
    a GMF (CDOP5n) for inversion.
    NOTE: wind and SAR datasets must be in a format supported by Nansat

    :param wind_data_src: /path/to/sar/data which can be opened with Nansat or Nansat object
    :param sar_data_src: /path/to/wind/data which can be opened with Nansat or Nansat object
    :param pixel_size_m: pixel size of the final product in meters
    """
    if isinstance(wind_source, (str, Path)):
        # Import and preprocess SAR data. 
        # NOTE 1: Only Sentinel-1 are supported at the moment
        # NOTE 2: To acquire better results do processing on full resolution and then 
        # resample inversed wind
        sar_data = preprocess_s1_data(sar_source, denoise_alg=denoise_alg, dst_px_size=pixel_size_m)
    else:
        sar_data = sar_source
    # Get arrays of SAR information required for the wind inversion: sigma0, incidence angle,
    # and look direction (required only for reprojecting wind direction)
    sigma0 = sar_data['sigma0_vv_denoised']
    incidence = sar_data['incidence_angle']
    look_dir = sar_data['look_direction']
    # Import wind model field, reproject, and extract wind speed and direction. 
    # NOTE 3: Only ERA5 is supported at the moment
    # NOTE 4: Wind speed is not used for the wind inversion
    aux_wind_spd, aux_wind_dir = preprocess_era5_data(wind_source, dst_geometry=sar_data)
    # Reproject wind direction to SAR antenna look direction
    aux_wind2sar_dir = wind2sar_direction(aux_wind_dir, look_dir)
    # Inverse wind speed from provided data
    sar_wind_spd = inverse_wind_spd(aux_wind2sar_dir, sigma0, incidence, 'CMOD5n')
    # Export product
    if export:
        _export2netcdf()

    return sar_wind_spd


def inverse_wind_spd(
        wind_dir: NDArray,
        sigma0: NDArray,
        incidence: NDArray,
        gmf_name: str = 'CMOD5n'
    ) -> NDArray:
    """
    Inverse wind speed from using CMOD5n GMF for given SAR sigma0, incidence,
    and polarizaion and collocated wind direction (e.g., from a model)
    
    :param wind_dir: 2D array, Wind direction with repect to SAR antenna look direction
        in degrees. Such as 0/180 deg. is wind blows towards/away from antenna
    :param sigma0: 
    :param incidence: 2D array, SAR incidence angle in degrees.
    :param gmf_name: Name of geophysical model function for the wind inversion
        NOTE: Only CMOD5n is supported at the moment
    :returns sar_wind_spd: Wind speed in m/s
    """
    if gmf_name == 'CMOD5n':
        sar_wind_spd = cmod5n_inverse(sigma0, wind_dir, incidence)
    else:
        raise ValueError
    return sar_wind_spd


def _export2netcdf():
    pass