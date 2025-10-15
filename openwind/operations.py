import logging
from collections.abc import Sequence
from datetime import datetime, timezone
from pathlib import Path
from typing import Union, Literal

import asf_search
import asf_search.constants.INTERNAL
import cartopy.crs as ccrs
import matplotlib.pyplot as plt
import netCDF4
import numpy as np
import xarray as xr

import openwind.gmf.cmod5n as cmod5n
import openwind.gmf.cmod7 as cmod7
import openwind.sar_data as sar_data
import openwind.utils as utils
import openwind.wind_data as wind_data


logger = logging.getLogger('openwind')
handler = logging.StreamHandler()
handler.setFormatter(logging.Formatter('%(asctime)s - %(name)s - %(levelname)s - %(message)s'))
logger.addHandler(handler)
logger.setLevel(logging.INFO)

asf_search.constants.INTERNAL.CMR_TIMEOUT = 60


def make_metadata(sar_source:sar_data.SARSource, wind_time, wind_model='ERA5', gmf='cmod5.n'):
    """Generate the metadata for the full product"""
    extent = sar_source.bounding_box
    time_start = sar_source.start_time
    time_end = sar_source.end_time

    return {
        'Conventions': 'CF-1.6',
        'institution': 'NERSC',
        'source': 'satellite remote sensing',
        'creation_date': datetime.now(timezone.utc).isoformat(),
        'northernmost_latitude': extent.north,
        'southernmost_latitude': extent.south,
        'westernmost_longitude': extent.west,
        'easternmost_longitude': extent.east,
        'BEAM_MODE': sar_source.properties['beamModeType'],
        'BEAM_SWATH': sar_source.properties['beamModeType'],
        'SWATH': sar_source.properties['beamModeType'],
        'data_center': '{"Bucket_Level0": "CONSORTIA/INSTITUTIONS", "Bucket_Level1": "", "Bucket_Level2": "", "Bucket_Level3": "", "Short_Name": "NERSC", "Long_Name": "Nansen Environmental and Remote Sensing Centre", "Data_Center_URL": "http://www.nersc.no/main/index2.php"}',
        'entry_title': 'Wind field from S1A_IW_GRDH_1SDV_20200404T044634_20200404T044659_031973_03B148_2646.nc',
        'instrument': '{"Category": "Earth Remote Sensing Instruments", "Class": "Active Remote Sensing", "Type": "Imaging Radars", "Subtype": "", "Short_Name": "SAR", "Long_Name": "Synthetic Aperture Radar"}',
        'ISO_topic_category': 'Imagery/Base Maps/Earth Cover',
        'keywords': "['Earth Science', 'Spectral/Engineering', 'RADAR', 'RADAR backscatter'], ['Earth Science', 'Spectral/Engineering', 'RADAR', 'RADAR imagery'], ['Earth Science', 'Spectral/Engineering', 'Microwave', 'Microwave Imagery'], ['EARTH SCIENCE', 'ATMOSPHERE', 'ATMOSPHERIC WINDS', 'SURFACE WINDS', 'U/V WIND COMPONENTS']",
        'keywords_vocabulary': 'GCMD Science Keywords',
        'MISSION_ID': sar_source.properties['fileID'].split('_')[0],
        'MODE': sar_source.properties['beamModeType'],
        'netcdf4_version_id': netCDF4.getlibversion().split()[0],
        'ORBIT_DIRECTION': sar_source.properties['flightDirection'],
        'ORBIT_NUMBER': sar_source.properties['orbit'],
        'platform': '{"Category": "Earth Observation Satellites", "Series_Entity": "SENTINEL-1", "Short_Name": "SENTINEL-1A", "Long_Name": "SENTINEL-1A"}',
        'polarisation': sar_source.properties['polarization'],
        'ProductTimelinessCategory': 'NTC',
        'PRODUCT_TYPE': sar_source.properties['processingLevel'],
        'SATELLITE_IDENTIFIER': sar_source.properties['platform'][:-1],
        'SENSOR_IDENTIFIER': sar_source.properties['sensor'],
        'summary': f'Near surface (10m) wind from Sentinel-1 C-band SAR GRD product and {wind_model} wind model, computed using the {gmf.upper()} algorithm (https://scatterometer.knmi.nl/)',
        'time_coverage_end': time_end.astimezone(timezone.utc).isoformat(),
        'time_coverage_start': time_start.astimezone(timezone.utc).isoformat(),
        'title': 'Near surface wind from Sentinel-1 GRD data',
        'winddir_time': wind_time.astimezone(timezone.utc).isoformat(),
        'WIND_DIRECTION_SOURCE': wind_model,
        'project': 'OpenWind',
    }


def make_full_dataset(sar_source: sar_data.SARSource, wind_source: wind_data.ERA5Source,
                      out_dir: Path,
                      gmf: Literal['cmod5.n', 'cmod7'] = 'cmod5.n', iterations: int = 10):
    """Generate the full dataset, which includes:
        - the original SAR data (sigma0, incidence angle, watermask)
        - the denoised SAR data (sigma0) if applicable
        - the original wind model data (speed+direction or u+v)
        - the wind speed generated from each sigma0 variable
    """
    out_file = Path(
        out_dir, f"{gmf.replace('.', '')}_{wind_source.model_name}_{sar_source.identifier}.nc")
    logger.info("Creating full dataset at %s", out_file)

    if out_file.exists():
        logger.info("Full dataset already exists at %s, skipping", out_file)
    else:
        s1_dataset = xr.open_dataset(sar_source.preprocessed_path, decode_coords='all')
        interp_wind_dataset = xr.open_dataset(wind_source.interpolated_path, decode_coords='all')

        wind_dir = interp_wind_dataset['dir'].to_numpy()

        if gmf == 'cmod5.n':
            gmf_inverse = cmod5n.cmod5n_inverse
        elif gmf == 'cmod7':
            gmf_inverse = cmod7.cmod7_inverse

        final_vars = {}
        for var_name in s1_dataset.variables:
            if var_name.startswith('sigma0'):
                wind_speed = gmf_inverse(
                    s1_dataset[var_name].to_numpy(),
                    wind_dir,
                    s1_dataset['incidence'].to_numpy(),
                    iterations=iterations
                )
                computed_u, computed_v = u_v_from_dir_speed(wind_dir, wind_speed)
                final_vars[var_name] = s1_dataset[var_name]
                final_vars[f"computed_wind_speed_from_{var_name}"] = (('row', 'col'), wind_speed)
                final_vars[f"computed_u_from_{var_name}"] = (('row', 'col'), computed_u)
                final_vars[f"computed_v_from_{var_name}"] = (('row', 'col'), computed_v)

        try:
            u10 = interp_wind_dataset['u10'].data
            v10 = interp_wind_dataset['v10'].data
        except KeyError:
            u10, v10 = u_v_from_dir_speed(interp_wind_dataset['dir'].data,
                                          interp_wind_dataset['speed'].data)
        final_vars["u10"] = (('row', 'col'), u10)
        final_vars["v10"] = (('row', 'col'), v10)

        full_dataset = xr.Dataset(
            attrs=make_metadata(
                sar_source,
                wind_source.get_time(),
                wind_source.model_name),
            coords={
                'lon': s1_dataset.coords['lon'],
                'lat': s1_dataset.coords['lat'],
            },
            data_vars={
                'watermask': s1_dataset['watermask'],
                'incidence': s1_dataset['incidence'],
                'model_wind_dir': interp_wind_dataset['dir'],
                'model_wind_speed': interp_wind_dataset['speed'],
                **final_vars,
            },
        )
        full_dataset.to_netcdf(out_file)
    return out_file


def u_v_from_dir_speed(dir, magnitude):
    """Calculates eastward and northward components from speed and direction
    """
    angle = np.deg2rad(90. - dir)
    return (-magnitude * np.cos(angle), -magnitude * np.sin(angle))


def plot_full_dataset(sar_source: sar_data.SARSource, wind_source,
                      full_ds_path: Path, out_dir: Path):
    """Makes a basic plot of the full dataset"""
    out_file = Path(out_dir, f'{full_ds_path.stem}.png')
    logger.info("Plotting full dataset. Writing to %s", out_file)
    if out_file.exists():
        logger.info("Plot already exists at %s, skipping", out_file)
    else:
        full_dataset = xr.open_dataset(full_ds_path, decode_coords='all')

        sigma0_variables = []
        computed_variables = []
        for variable in full_dataset.variables:
            if variable.startswith('sigma0'):
                sigma0_variables.append(variable)
            elif variable.startswith('computed_'):
                computed_variables.append(variable)
        sigma0_variables.sort()
        computed_variables.sort()
        base_variables = set(('u10', 'v10', 'model_wind_dir', 'model_wind_speed'))

        # normalize s1
        watermask = full_dataset['watermask'].to_numpy()
        mask = None
        for sigma0_var in sigma0_variables:
            denoised = full_dataset[sigma0_var].to_masked_array(copy=False)
            mask = np.logical_or(denoised.mask, watermask != 1.)
            flat = denoised[~mask].flatten()
            low, high = np.percentile(flat, (5, 95))
            np.clip(denoised, low, high, denoised, casting='unsafe')
            #normalize in the interval [0,100]
            denoised[~mask] = ((denoised[~mask] - low) / (high - low)) * 100.
            denoised[mask] = np.nan

        # mask other variables
        for variable in (*base_variables, *computed_variables):
            full_dataset[variable].to_masked_array(copy=False)[mask] = np.nan

        nr_sigma0 = len(sigma0_variables)

        ncols = nr_sigma0 if nr_sigma0 >= 2 else 2
        nrows = 2 + nr_sigma0

        fig, axs = plt.subplots(
            subplot_kw={'projection': ccrs.epsg(3857)},
            ncols=ncols,
            nrows=nrows,
            squeeze=False,
            figsize=(2.*ncols, 2.*nrows),
        )

        for row in axs:
            for ax in row:
                ax.set_extent(sar_source.bounding_box)
                ax.coastlines()
                ax.set_anchor('W')
                ax.set_xlabel('')
                ax.set_ylabel('')

        for i, var in enumerate(sigma0_variables):
            full_dataset[var].plot(
                ax=axs[0, i],
                x='lon', y='lat', cmap='gray', transform=ccrs.PlateCarree(),
                add_labels=False, add_colorbar=False)
            axs[0, i].set_title(var, fontsize='small')

            full_dataset.plot.streamplot(
                ax=axs[2+i, 0], zorder=3,
                x='lon', y='lat', u=f'computed_u_from_{var}', v=f'computed_v_from_{var}',
                transform=ccrs.PlateCarree(),
                linewidth=.5, arrowsize=.5, density=3, hue=f'computed_wind_speed_from_{var}',
                cbar_kwargs={'ax': axs[2+i, 0], 'label': '', 'shrink': .7, 'extend': 'both', 'location': 'right'})
            axs[2+i, 0].set_title(f"Wind from {var}", fontsize='small')

            full_dataset[f'computed_wind_speed_from_{var}'].plot(
                ax=axs[2+i, 1],
                x='lon', y='lat', transform=ccrs.PlateCarree(),
                add_labels=False, #add_colorbar=False,
                cbar_kwargs={'ax': axs[2+i, 1], 'label': '', 'shrink': .7, 'extend': 'both', 'location': 'right'})
            axs[2+i, 1].set_title(f"Wind speed from {var}", fontsize='small')

        full_dataset.plot.streamplot(
            ax=axs[1, 0], zorder=3,
            x='lon', y='lat', u='u10', v='v10', transform=ccrs.PlateCarree(),
            linewidth=.5, arrowsize=.5, density=3, hue='model_wind_speed',
            cbar_kwargs={'ax': axs[1, 0], 'label': '', 'shrink': .7, 'extend': 'both', 'location': 'right'})
        axs[1, 0].set_title(f'{wind_source.model_name} wind', fontsize='small')

        full_dataset['model_wind_speed'].plot(
            ax=axs[1, 1],
            x='lon', y='lat', transform=ccrs.PlateCarree(),
            add_labels=False, #add_colorbar=False,
            cbar_kwargs={'ax': axs[1, 1], 'label': '', 'shrink': .7, 'extend': 'both', 'location': 'right'})
        axs[1, 1].set_title(f"{wind_source.model_name} wind speed", fontsize='small')

        fig.tight_layout()
        fig.savefig(out_file, dpi=300)
    return out_file


def check_sar_input(
        extent: utils.Extent = None,
        time_start: datetime = None,
        time_end: datetime = None,
        s1_identifiers: Sequence[str] = None,
        sar_files: Sequence[Union[str, Path]] = None):
    """Check that input parameters give enough information to look for
    SAR sources
    """
    if sar_files is None and s1_identifiers is None:
        for param in (extent, time_start, time_end):
            if param is None:
                raise ValueError(
                    "Either sar_files, s1_identifiers or (extent, time_start, time_end) "
                    "need to be provided")
            else:
                sar_input_params = f"extent={extent}, time_start={time_start}, time_end={time_end}"
    else:
        sar_input_params = f"sar_files={sar_files}"
    return sar_input_params


def create_sar_sources(
        extent: utils.Extent = None,
        time_start: datetime = None,
        time_end: datetime = None,
        s1_identifiers: Sequence[str] = None,
        sar_files: Sequence[Union[str, Path]] = None,
        sar_source_class: Union[type[sar_data.Sentinel1Source],
                                type[sar_data.EnvisatASARSource]] = sar_data.Sentinel1Source):
    """Create SAR sources from input parameters"""
    sar_sources = None
    if s1_identifiers:
        sar_sources = sar_data.Sentinel1Source.from_asf(
            query={
                'granule_list': s1_identifiers,
                'processingLevel': asf_search.PRODUCT_TYPE.GRD_HD,
            }
        )
    elif sar_source_class is not None and sar_files is not None:
        sar_sources = [sar_source_class.from_path(data_path=p) for p in sar_files]
    else:
        sar_sources = sar_data.Sentinel1Source.from_asf(extent, time_start, time_end)
    if sar_sources is None:
        raise ValueError("Could not determine the source of SAR data")
    elif not sar_sources:
        raise RuntimeError(
            f"Could not find any SAR sources matching the provided parameters ({sar_input_params})")
    return sar_sources


def generate_product(
        extent: utils.Extent = None,
        time_start: datetime = None,
        time_end: datetime = None,
        s1_identifiers: Sequence[str] = None,
        sar_files: Sequence[Union[str, Path]] = None,
        sar_source_class: Union[type[sar_data.Sentinel1Source],
                                type[sar_data.EnvisatASARSource]] = sar_data.Sentinel1Source,

        wind_source_class: Union[type[wind_data.ERA5Source],
                                 type[wind_data.Sentinel1OCNSource]] = wind_data.ERA5Source,
        gmf: Literal['cmod5.n', 'cmod7'] = 'cmod5.n',
        iterations: int = 10,
        plot: bool = False,
        workdir: Union[str, Path] = Path('.'),
        input_dir: Union[str, Path] = None,
        output_dir: Union[str, Path] = None,
        denoised_dir: Union[str, Path] = None,
        wind_folder: Union[str, Path] = None,
        plot_dir: Union[str, Path] = None):
    """Generate a full data product for the specified time and space extents"""
    sar_input_params = check_sar_input(extent, time_start, time_end, s1_identifiers, sar_files)

    # Make sure necessary directories exist
    workdir = Path(workdir)
    input_dir = input_dir or (workdir / 'input')
    output_dir = output_dir or (workdir / 'output')
    denoised_dir = denoised_dir or (workdir / 'denoised')
    wind_folder = wind_folder or (workdir / 'wind')
    plot_dir = plot_dir or (workdir / 'plots')
    for folder in (input_dir, output_dir, denoised_dir, wind_folder, plot_dir):
        logger.debug("Making sure directory exists: %s", folder)
        folder.mkdir(parents=True, exist_ok=True)

    sar_sources = create_sar_sources(
        extent, time_start, time_end,
        s1_identifiers,
        sar_files, sar_source_class)

    for sar_source in sar_sources:
        logger.info("Processing %s", sar_source.identifier)

        # download SAR data
        sar_source.download(input_dir, unzip=True)
        # denoise and/or resize SAR data to a 500m grid
        sar_source.preprocess(denoised_dir, polarizations=('VV',))

        # download wind data
        wind_source = wind_source_class(sar_source)
        try:
            wind_source.download(wind_folder)
        except RuntimeError:
            logger.warning("No wind data found matching %s", sar_source.identifier)
            continue

        # interpolate wind data on the SAR grid
        wind_source.interpolate_on_sar_grid(wind_folder)

        # create the full dataset
        full_dataset_path = make_full_dataset(
            sar_source, wind_source, output_dir, gmf, iterations=iterations)

        # plot the full dataset
        if plot:
            plot_full_dataset(sar_source, wind_source, full_dataset_path, plot_dir)

    logger.info('Done processing sar sources matching %s', sar_input_params)
