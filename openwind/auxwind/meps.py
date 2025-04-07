from netCDF4 import Dataset, num2date
from datetime import datetime, timedelta
import numpy as np
import os
from nansat import Nansat, Domain
import pandas as pd
import xarray as xr

class MEPSFilename:
    
    FILENAME_PATTERN = 'meps_lagged_6_h_subset_2_5km_%(year)d%(month).2d%(day).2dT%(hour).2dZ.nc'
    FILE_SRC_PATTERN = 'https://thredds.met.no/thredds/dodsC/meps25epsarchive/%(year)d/%(month).2d/%(day).2d'

    def __init__(self, target_datetime):
        self.target_datetime = target_datetime
        self.meps_datetime = self.estimate_meps_datetime()
        self.meps_url = self.assemble_filename()
        
    def estimate_meps_datetime(self):
        """The model runs at 00 06 12 18 UTC therefore we need to 
        find last hour before Sentinel-1 acquisition"""
        meps_hours = [0, 6, 12, 18]
        year = self.target_datetime.year
        month = self.target_datetime.month
        day = self.target_datetime.day
        # Generate a list of datetime candidates for the day of SAR acquisition
        meps_datetimes = np.array([datetime(year, month, day, hour, 0, 0) 
                                   for hour in meps_hours])
        # Find time difference between each candidate and target_time 
        dt = meps_datetimes - self.target_datetime
        # Find closest timedelta before the target_time
        closest_timedelta_before = dt[[dti.days < 0 for dti in dt]].max()
        meps_time = self.target_datetime + closest_timedelta_before
        return meps_time
    
    def assemble_filename(self):
        meps_src_raw = os.path.join(self.FILE_SRC_PATTERN, self.FILENAME_PATTERN)
        meps_src = meps_src_raw % {'year': self.meps_datetime.year, 
                                   'month': self.meps_datetime.month, 
                                   'day': self.meps_datetime.day, 
                                   'hour': self.meps_datetime.hour}
        return meps_src


def lcc_grid_convergence_angle(
        longitude: float, 
        central_meridian: float = -25., 
        central_latitude: float = 77.5) -> float:
    """ 
    Calculate grid convergence angle for the Lambert Conic Conformal (LCC)
    projection

    :param longitude:        Longitude of the grid point from the LCC domain
    :param central_meridian: Central meridial of the LCC grid projection
    :param central_latitude: Standard parallel of the LCC grid projection
    :returns azimuth:        Deviation of the LCC domain from the true north
    """
    azimuth = np.sin(np.deg2rad(central_latitude)) * (longitude - central_meridian)
    return azimuth


def calculate_wind_magnitude(u, v):
    """
    Calculate wind speed from the components
    
    :param u: Eastward wind component in m/s
    :param v: Northward wind component in m/s
    :returns float: Return wind speed
    """
    return np.hypot(u, v)


def calculate_wind_direction(u, v):
    """
    Estimate wind direction from components. Meteo convention
    
    :param u: Eastward wind component in m/s
    :param v: Northward wind component in m/s
    :returns float: Return wind direction in degrees 
    """
    direction = (-np.rad2deg(np.arctan2(v, u)) - 90)
    direction[direction < 0] += 360
    return direction


def get_uv_wind_components(magnitude, direction):
    """
    Decompose wind vector for the eastwards and northwards components
    :param magnitude: Wind speed in m/s
    :param direction: Wind direction in degrees
    :return u, v: eastward and northward wind components
    """
    u = magnitude * np.cos(np.deg2rad(-direction - 90))
    v = magnitude * np.sin(np.deg2rad(-direction - 90))
    return u, v


def process_wind(dataset):
    # Process wind parameters from the x and y components
    # 1. xw wind speed and direction from the components
    dataset['wind_speed'] = xr.Variable(
        ('ensemble_member', 'y', 'x'),
        calculate_wind_magnitude(dataset['x_wind_10m'].data, dataset['y_wind_10m'].data))
    
    wind_direction = calculate_wind_direction(dataset['x_wind_10m'].data, dataset['y_wind_10m'].data)
    # 2. compensate on grid angle deviation from true north
    # Calculate grid convergence correction
    grid_convergence_angle = lcc_grid_convergence_angle(dataset['longitude'].data)
    dataset['wind_direction'] = xr.Variable(
        ('ensemble_member', 'y', 'x'),
        (wind_direction + grid_convergence_angle) % 360)
    # 3. retrieve u and v component
    u, v = get_uv_wind_components(dataset['wind_speed'].data, dataset['wind_direction'].data)
    dataset['u'] = xr.Variable(('ensemble_member', 'y', 'x'), u)
    dataset['v'] = xr.Variable(('ensemble_member', 'y', 'x'), v)
    return dataset


class MEPSDataOpenDAP:
    
    EXTENT_PATTERN = '-te %(min_x)f %(min_y)f %(max_x)f %(max_y)f -ts %(x_size)d %(y_size)d'
    
    def __init__(self, src):
        self.src = src
        self.netcdf_dataset = Dataset(src)
        # Get model spatial domain
        self.domain = self.generate_domain()
        # Calculate grid convergence angle
        self.grid_convergence_angle = self.calculate_grid_convergence_angle()

    def generate_domain(self):
        domain_srs = self.netcdf_dataset['projection_lambert'].proj4
        domain_ext = self.get_domain_extent()
        d = Domain(domain_srs, domain_ext)
        return d
    
    def get_domain_extent(self):
        xs = self.netcdf_dataset['x'][:]
        ys = self.netcdf_dataset['y'][:]
        # Fill pattern with the data from the dataset
        domain_extent = self.EXTENT_PATTERN % {'min_x': xs.min(), 'max_x': xs.max(),
                                               'min_y': ys.min(), 'max_y': ys.max(), 
                                               'x_size': xs.size, 'y_size': ys.size}
        return domain_extent
    
    def calculate_grid_convergence_angle(self):
        meps_grid_convergence_angle = grid.lcc_grid_convergence_angle(
            longitude=self.netcdf_dataset['longitude'][:].data,
            central_meridian=self.netcdf_dataset['projection_lambert'].longitude_of_central_meridian,
            central_latitude=self.netcdf_dataset['projection_lambert'].latitude_of_projection_origin
        )
        return meps_grid_convergence_angle

    
    def find_closest_forecast_time(self, target_time):
        # Find closest model forecast to the sar acquisition
        meps_forecasts = num2date(self.netcdf_dataset['time'][:], 
                                  self.netcdf_dataset['time'].units)
        dt = meps_forecasts - target_time
        closest_after_sar = dt[[dti.days >= 0 for dti in dt]][0].seconds
        closest_before_sar = dt[[dti.days < 0 for dti in dt]][-1].seconds 
        # compensate for the -1 day
        closest_before_sar = 24 * 60 * 60 - closest_before_sar
        # Chhose the closes forecast and retrieve the full datetime for it
        if closest_before_sar <= closest_after_sar:
            closest_forecast_time = target_time - timedelta(seconds=closest_before_sar)
        else:
            closest_forecast_time = target_time + timedelta(seconds=closest_after_sar)
            
        closest_forecats_time_id = np.where(meps_forecasts == closest_forecast_time)[0][0]
        return closest_forecast_time, closest_forecats_time_id
    
    
    def export_wind_forecast(self, target_time, export=True, dst=''):
        # Find closest in time forecast in the file
        forecast_time, forecast_id = self.find_closest_forecast_time(target_time)
        # Czlculate wind field parameters (true North)
        wind_speed_m, wind_direction_m, u_m, v_m = self.process_wind(forecast_id, ensemble_mean=True)
        wind_speed_d, wind_direction_d, u_d, v_d = self.process_wind(forecast_id, deterministic=True)
        if export:
            meps_data = Nansat.from_domain(self.domain)
            # Set metadata from the model time
            meps_metadata = {att:self.netcdf_dataset.getncattr(att) for 
                             att in self.netcdf_dataset.ncattrs()}
            meps_data.set_metadata(meps_metadata)
            # Add grid convergence angle 

            # Set time coverage start/end
            meps_data.set_metadata(key='time_coverage_start', value=forecast_time)
            meps_data.set_metadata(key='time_coverage_end', value=forecast_time)

            # Add bands
            # Add ensemble mean estimates
            meps_data.add_band(wind_speed_m, parameters={'name': 'wind_speed_m', 'units': 'm s-1',
                                                         'long_name': 'ensemble_mean_wind_speed_at_10m_height'})
            meps_data.add_band(wind_direction_m, parameters={'name': 'wind_direction_m', 'units': 'degree',
                                                             'long_name': 'ensemble_mean_wind_direction_at_10m_height', 
                                                             'standard_name': 'wind_from_direction'})
            meps_data.add_band(u_m, parameters={'name': 'u_m', 'units': 'm s-1', 'standard_name': 'eastward_wind', 
                                                'long_name': 'ensemble_mean_eastward_wind'})
            meps_data.add_band(v_m, parameters={'name': 'v_m', 'units': 'm s-1','standard_name': 'northward_wind',
                                                'long_name': 'ensemble_mean_northward_wind'})
            # Add determenistic estimates
            meps_data.add_band(wind_speed_d, parameters={
                'name': 'wind_speed_d', 'units': 'm s-1', 
                'long_name': 'deterministic_wind_speed_at_10m_height'})
            meps_data.add_band(wind_direction_d, parameters={
                'name': 'wind_direction_d', 'units': 'degree', 'standard_name': 'wind_from_direction',
                'long_name': 'deterministic_wind_direction_at_10m_height'})
            meps_data.add_band(u_d, parameters={
                'name': 'u_d', 'units': 'm s-1', 'standard_name': 'eastward_wind', 
                'long_name': 'deterministic_eastward_wind'})
            meps_data.add_band(v_d, parameters={
                'name': 'v_d', 'units': 'm s-1','standard_name': 'northward_wind',
                'long_name': 'deterministic_northward_wind'})

            if len(dst) == 0:
                dst = './meps_subset_2_5km_' + forecast_time.isoformat() + '.nc'
                dst = dst.replace('-', '').replace(':', '')
            
            # if reproject2dom:
            #     meps_data.reproject(reproject2dom)

            # meps_data.export(dst)
            return meps_data
