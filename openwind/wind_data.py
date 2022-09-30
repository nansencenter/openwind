# -------------------------------------------------------------------------------
# Name:		wind_data.py
# Purpose:  Collection of utility functions to access and preprocess auxiliary 
# wind information from numerical models such as ERA5 and NCEP GFS 
#
# Author:   Artem Moiseev
# Modified:	
#
# Created:	24.09.2022
# Last modified: 
# Copyright: (c) NERSC
# License:      
# -------------------------------------------------------------------------------
import numpy as np
import cdsapi
from nansat import Nansat
from pathlib import Path
from typing import Union, Optional, Tuple
from numpy.typing import NDArray
from datetime import datetime


def direction_from(u, v):
    """
    Estimate direction from components.
    
    NOTE: Meteorological convention (i.e. direction_from)
    :param u: Eastward component in m/s
    :param v: Northward component in m/s
    :returns direction: Direction in degrees 
    """
    direction = (-np.rad2deg(np.arctan2(v, u)) - 90) % 360
    return direction


def magnitude(u, v):
    """
    Calculate vector magnitude from the components
    
    :param u: Eastward component in m/s
    :param v: Northward component in m/s
    :returns: Magnitude in m/s
    """
    return np.hypot(u, v)


def wind2sar_direction(
        wind_dir: Union[int, float, NDArray],
        look_dir: Union[int, float, NDArray]
    ) -> Union(NDArray, float):
    """
    Reproject wind direction from geographic (where 0 deg. is northward wind) to SAR
    antenna look direction (where 0 deg. is wind towards SAR)

    :param wind_dir: Wind direction in geographic system (0 deg. is north)
    :param look_dir: SAR antenna look direction
    :returns: Wind direction in deg where 0/180 deg is wind toward/away from antenna
    """
    return np.mod(wind_dir - look_dir, 360)


def fetch_era5_data(
        timestamp: datetime,
        central_lat: Union(int, float),
        central_lon: Union(int, float),
        dst: Union(str, Path),
        pad: Union(int, float) = 5
    ) -> Path:
    """
    Download ERA5 reanalysis wind field at 10 m (u and v components) from the Copernicus 
    Climate Change Service (https://cds.climate.copernicus.eu/cdsapp#!/dataset/reanalysis
    -era5-singale-levels?tab=overview) using CDS API
    
    :param timestamp: datetime object
    :param central_lat: central latitude of the area
    :param central_lon: central longitude of the area
    :param pad: distance from the center of the search area
    :returns dst_uri: Path to the downloaded NetCDF file
    """
    # Configure time and domain for the cds api query
    cds_query = {
        'product_type': 'reanalysis',
        'format': 'netcdf',
        'variable': ['10m_u_component_of_wind', '10m_v_component_of_wind'],
        'month': f'{timestamp.month:02d}',
        'year': f'{timestamp.year}',
        'day': f'{timestamp.day:02d}',
        'time': f'{timestamp:%H}:00',
        'area': [central_lat - pad, central_lon - pad, 
                 central_lat + pad, central_lon + pad]}
    # Connect to the client
    cds_client = cdsapi.Client()
    # Create dst uri
    dst_uri = Path(dst) / f'ERA5_{timestamp:%Y%m%dT%H}00.nc'
    print(f'>> Downloading {dst_uri}')
    # Download the data
    cds_client.retrieve('reanalysis-era5-single-levels', cds_query, dst_uri)
    return dst_uri


def preprocess_wind_data(
        source: Union[str, Path, Nansat], 
        dst_geometry: Optional[Nansat] = None
    ) -> Tuple[NDArray, NDArray]:
    """
    Read wind dataset and extract wind speed and direction.
    NOTE: Only ERA5 is supported at the moment
    
    :param uri: /path/to/wind/file that can be opened with nansat or Nansat object
    :param dst_geometry: destination grid. Bilinear resampling is applied.
        geometry can me S1 frame opened in Nansat or nansat Domain generated from lon lat grids
    :returns wind_spd, wind_dir: numpy arrays with wind speed in m/s and direction in deg
        (direction from)
    """
    # If source is path to the file then read the data using Nansat
    if isinstance(source, (str, Path)):
        print(f'>> Reading {source}')
        # Read era5 data using Nansat
        wind_data = Nansat(str(source))
    # If source is nansat type object then use it for the processing directly 
    else:
        wind_data = source
    # If dst geometry provided then reproject era5 data to the dst geometry
    if dst_geometry is not None:
        # Resample using bilinear interpolation
        wind_data.reproject(dst_geometry, resample_alg=1)
    # Calculate wind speed and direction from u and v components provided in model
    wind_spd = magnitude(wind_data['u10'], wind_data['v10'])
    wind_dir = direction_from(wind_data['u10'], wind_data['v10'])
    return wind_spd, wind_dir
