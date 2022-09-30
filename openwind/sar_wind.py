from numpy.typing import NDArray
from typing import Union, Optional, Tuple
from nansat import Nansat
from openwind.sar_data import preprocess_sar_data
from openwind.wind_data import preprocess_era5_data, wind2sar_direction
from openwind.gmf.cmod5n import cmod5n_inverse
from pathlib import Path
import numpy as np


def derive_sar_wind(
        wind_source: Union[str, Path, Nansat],
        sar_source: Union[str, Path, Nansat],
        pixel_size_m: Union[int, float] = None,
        denoise_alg: Optional[str] = None,
        export_dst: Optional[Path] = None
    ) -> Tuple[NDArray, NDArray]:
    """
    Derive wind product from collocated SAR observations and model field using
    a GMF (CDOP5n) for inversion.
    NOTE: wind and SAR datasets must be in a format supported by Nansat

    :param wind_data_src: /path/to/sar/data which can be opened with Nansat or Nansat object
    :param sar_data_src: /path/to/wind/data which can be opened with Nansat or Nansat object
    :param pixel_size_m: pixel size of the final product in meters
    :param denoise_alg: Denoising algorithm name
        NOTE: Only applicable to Sentinel-1 data
    :param export_dst: export wind product to a netcdf
    :returns sar_wind_spd, aux_wind_dir: SAR derived wind speed in m/s and aux wind direction 

    """
    if isinstance(sar_source, (str, Path)):
        # Ensure that the path to the source file is Path object
        sar_source = Path(sar_source)
        # Import and preprocess SAR data. 
        # NOTE: To acquire better results do processing on full resolution and then 
        # resample inversed wind
        sar_data = preprocess_sar_data(sar_source, denoise_alg=denoise_alg, dst_px_size=pixel_size_m)
    else:   
        sar_data = sar_source

    if Path(sar_data.filename).name.startswith('S1') and denoise_alg is not None:
        sigma0_bandname = 'sigma0_vv_denoised'
    else:
        sigma0_bandname = 'sigma0_VV'
    # Get arrays of SAR information required for the wind inversion: sigma0, incidence angle,
    # and look direction (required only for reprojecting wind direction)
    sigma0 = sar_data[sigma0_bandname]
    incidence = sar_data['incidence_angle']
    look_dir = sar_data['look_direction']
    # Import wind model field, reproject, and extract wind speed and direction. 
    # NOTE: Only ERA5 is supported at the moment
    # NOTE: Wind speed is not used for the wind inversion
    aux_wind_spd, aux_wind_dir = preprocess_era5_data(wind_source, dst_geometry=sar_data)
    # Reproject wind direction to SAR antenna look direction
    aux_wind2sar_dir = wind2sar_direction(aux_wind_dir, look_dir)
    # Inverse wind speed from provided data
    sar_wind_spd = inverse_wind_spd(aux_wind2sar_dir, sigma0, incidence, 'CMOD5n')
    # Export product
    if export_dst is not None:
        _export2netcdf()

    return sar_wind_spd, aux_wind_dir


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


def _export2netcdf(export_dst):
    pass