import logging
import multiprocessing
import shutil
from collections.abc import Sequence
from concurrent.futures import ThreadPoolExecutor, as_completed
from datetime import datetime, timezone
from pathlib import Path
from typing import Union, Literal, Generator

import asf_search
import asf_search.constants.INTERNAL
import cartopy.crs as ccrs
import cartopy.feature.download.__main__
import matplotlib.pyplot as plt
import netCDF4
import numpy as np
import requests
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


class Done:
    """Put in a Queue to signal the end of processing"""


def make_metadata(sar_source:sar_data.SARSource,
                  wind_time, wind_model='ERA5',
                  gmf='cmod5.n',
                  pixel_size: int = 500):
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
        'data_center': "Nansen Environmental and Remote Sensing Centre",
        'title': f'Wind field from {sar_source.identifier}',
        'instrument': 'SAR (Synthetic Aperture Radar)',
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
        'summary': (
            f'Near surface (10m) wind from SAR ({sar_source.platform}) and {wind_model} wind model,'
            ' computed using the {gmf.upper()} GMF (https://scatterometer.knmi.nl/). This product '
            'was generated using Openwind (https://github.com/nansencenter/openwind/)'),
        'time_coverage_end': time_end.astimezone(timezone.utc).isoformat(),
        'time_coverage_start': time_start.astimezone(timezone.utc).isoformat(),
        'title': 'Near surface wind from Sentinel-1 GRD data',
        'winddir_time': wind_time.astimezone(timezone.utc).isoformat(),
        'WIND_DIRECTION_SOURCE': wind_model,
        'project': 'C3-eKerala',
        'pixel_size': f"{pixel_size}m"
    }


def make_full_dataset(sar_source: sar_data.SARSource, wind_source: wind_data.ERA5Source,
                      out_dir: Path,
                      full_ds_file_name: str = None,
                      gmf: Literal['cmod5.n', 'cmod7'] = 'cmod5.n', iterations: int = 10):
    """Generate the full dataset, which includes:
        - the original SAR data (sigma0, incidence angle, watermask)
        - the denoised SAR data (sigma0) if applicable
        - the original wind model data (speed+direction or u+v)
        - the wind speed generated from each sigma0 variable
    """
    if not full_ds_file_name:
        full_ds_file_name = (
            f"{gmf.replace('.', '')}_{wind_source.model_name}_{sar_source.identifier}.nc")
    out_file = Path(out_dir, full_ds_file_name)
    logger.info("Creating full dataset at %s", out_file)

    if out_file.exists():
        logger.info("Full dataset already exists at %s, skipping", out_file)
    else:
        with xr.open_dataset(sar_source.preprocessed_path, decode_coords='all') as s1_dataset, \
             xr.open_dataset(wind_source.interpolated_path, decode_coords='all'
                             ) as interp_wind_dataset:

            wind_sar_dir = interp_wind_dataset['sar_dir'].to_masked_array(copy=False)
            dimensions = ('time', 'row', 'col')

            if gmf == 'cmod5.n':
                gmf_inverse = cmod5n.cmod5n_inverse
            elif gmf == 'cmod7':
                gmf_inverse = cmod7.cmod7_inverse

            sar = s1_dataset['sigma0'].to_numpy()
            # do not generate wind speed where the wind direction
            # is not available
            sar[wind_sar_dir.mask] = np.nan
            wind_speed = gmf_inverse(
                sar,
                wind_sar_dir,
                s1_dataset['angle_of_incidence'].to_numpy(),
                iterations=iterations
            )

            full_dataset = xr.Dataset(
                attrs=make_metadata(
                    sar_source,
                    wind_source.get_time(),
                    wind_source.model_name,
                    s1_dataset.attrs['pixel_size']),
                coords={
                    'longitude': (dimensions, [s1_dataset.coords['longitude'].data]),
                    'latitude': (dimensions, [s1_dataset.coords['latitude'].data]),
                    'time': [np.datetime64(sar_source.start_time.isoformat())],
                },
                data_vars={
                    'sea_binary_mask': (dimensions, [s1_dataset['sea_binary_mask'].data]),
                    'angle_of_incidence': (dimensions, [s1_dataset['angle_of_incidence'].data]),
                    'direction_of_radial_vector_away_from_instrument': (
                        dimensions, [s1_dataset['look_direction'].data]),
                    'wind_from_direction_model': (dimensions, [interp_wind_dataset['dir'].data]),
                    'wind_speed_model': (dimensions, [interp_wind_dataset['speed'].data]),
                    'sigma0': (dimensions, [s1_dataset['sigma0'].data], s1_dataset['sigma0'].attrs),
                    'wind_speed_sar': (dimensions, [wind_speed]),
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
        with xr.open_dataset(full_ds_path, decode_coords='all') as full_dataset, \
             xr.open_dataset(wind_source.interpolated_path, decode_coords='all'
                             ) as interp_wind_dataset:
            full_ds_time0 = full_dataset.isel(time=0)

            # normalize s1
            watermask = full_ds_time0['sea_binary_mask'].to_numpy()
            mask = None
            sigma0_var = 'sigma0'
            sigma0_data = full_ds_time0[sigma0_var].to_masked_array(copy=False)
            mask = np.logical_or(sigma0_data.mask, watermask != 1.)
            flat = sigma0_data[~mask].flatten()
            low, high = np.percentile(flat, (5, 95))
            np.clip(sigma0_data, low, high, sigma0_data, casting='unsafe')
            #normalize in the interval [0,100]
            sigma0_data[~mask] = ((sigma0_data[~mask] - low) / (high - low)) * 100.
            sigma0_data[mask] = np.nan

            # mask other variables
            min_wind = 0.
            max_wind = 20.

            fig, axs = plt.subplots(
                subplot_kw={'projection': ccrs.epsg(3857)},
                ncols=2,
                nrows=3,
                squeeze=False,
                figsize=(4., 6.),
            )

            for row in axs:
                for ax in row:
                    ax.set_extent(sar_source.bounding_box)
                    ax.coastlines()
                    ax.set_anchor('W')
                    ax.set_xlabel('')
                    ax.set_ylabel('')

            # SAR data
            full_ds_time0[sigma0_var].plot(
                ax=axs[0, 0],
                x='longitude', y='latitude', cmap='gray', transform=ccrs.PlateCarree(),
                add_labels=False, add_colorbar=False)
            axs[0, 0].set_title(sigma0_var, fontsize='small')

            # wind streamlines: direction from model, speed from SAR
            u, v = u_v_from_dir_speed(
                full_ds_time0['wind_from_direction_model'].to_masked_array(copy=True),
                full_ds_time0['wind_speed_sar'].to_masked_array(copy=True))
            u[mask] = np.nan
            v[mask] = np.nan
            temp_ds = xr.Dataset(
                coords=full_ds_time0.coords,
                data_vars={
                    'u': (full_ds_time0.dims, u),
                    'v': (full_ds_time0.dims, v),
                    'wind_speed_sar': full_ds_time0['wind_speed_sar']
                }
            )
            temp_ds.plot.streamplot(
                ax=axs[1, 0], zorder=3,
                x='longitude', y='latitude', u='u', v='v',
                vmin=min_wind, vmax=max_wind,
                transform=ccrs.PlateCarree(),
                linewidth=.5, arrowsize=.5, density=3, hue=f'wind_speed_sar',
                cbar_kwargs={'ax': axs[1, 0], 'label': '', 'shrink': .5,
                                'extend': 'both', 'location': 'right'})
            axs[1, 0].set_title(f"Wind from SAR", fontsize='small')

            # wind speed from SAR
            full_ds_time0['wind_speed_sar'].to_masked_array(copy=False)[mask] = np.nan
            full_ds_time0['wind_speed_sar'].plot(
                ax=axs[1, 1],
                x='longitude', y='latitude', transform=ccrs.PlateCarree(),
                vmin=min_wind, vmax=max_wind,
                add_labels=False, #add_colorbar=False,
                cbar_kwargs={'ax': axs[1, 1], 'label': '', 'shrink': .5,
                                'extend': 'both', 'location': 'right'})
            axs[1, 1].set_title("Wind speed from SAR", fontsize='small')

            # wind streamlines from model
            interp_wind_dataset.plot.streamplot(
                ax=axs[2, 0], zorder=3,
                x='longitude', y='latitude', u='u10', v='v10', transform=ccrs.PlateCarree(),
                vmin=min_wind, vmax=max_wind,
                linewidth=.5, arrowsize=.5, density=3, hue='speed',
                cbar_kwargs={'ax': axs[2, 0], 'label': '', 'shrink': .5,
                             'extend': 'both', 'location': 'right'})
            axs[2, 0].set_title(f'{wind_source.model_name} wind', fontsize='small')

            # wind speed from model
            full_ds_time0['wind_speed_model'].to_masked_array(copy=False)[mask] = np.nan
            full_ds_time0['wind_speed_model'].plot(
                ax=axs[2, 1],
                x='longitude', y='latitude', transform=ccrs.PlateCarree(),
                vmin=min_wind, vmax=max_wind,
                add_labels=False, #add_colorbar=False,
                cbar_kwargs={'ax': axs[2, 1], 'label': '', 'shrink': .5, 'extend': 'both', 'location': 'right'})
            axs[2, 1].set_title(f"{wind_source.model_name} wind speed", fontsize='small')

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


def make_product_file_name(sar_source: sar_data.SARSource, product_version: str, file_version: str):
    """Build the name for a final product file"""
    return (
        f"{sar_source.platform}_{sar_source.start_time.strftime('%Y%m%d%H%M%S')}"
        f"_SARWIND_v{product_version}_fv{file_version}.nc")


def check_final_file(out_dir: Union[str, Path],
                     sar_source: sar_data.SARSource,
                     product_version: str, file_version: str,
                     remove_invalid: bool = True):
    """Checks if the product file has already been generated and
    removes it if invalid
    """
    full_path = Path(out_dir, make_product_file_name(sar_source, product_version, file_version))
    result = False
    if full_path.is_file():
        try:
            netCDF4.Dataset(full_path)
        except OSError:
            if remove_invalid:
                logger.warning("Removing invalid product file: %s", full_path)
                full_path.unlink()
            else:
                logger.error(
                    "Won't process SAR source, invalid product file already present: %s", full_path)
        else:
            logger.info("Product file already exists: %s", full_path)
            result = full_path
    return result


def create_sar_sources(
        extent: utils.Extent = None,
        time_start: datetime = None,
        time_end: datetime = None,
        s1_identifiers: Sequence[str] = None,
        sar_files: Sequence[Union[str, Path]] = None,
        sar_source_class: Union[type[sar_data.Sentinel1Source],
                                type[sar_data.EnvisatASARSource]] = sar_data.Sentinel1Source,
        sar_input_params: str = None) -> Generator:
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
        sar_sources = (sar_source_class.from_path(data_path=p) for p in sar_files)
    else:
        sar_sources = sar_data.Sentinel1Source.from_asf(extent, time_start, time_end)
    if sar_sources is None:
        raise ValueError("Could not determine the source of SAR data")
    elif not sar_sources:
        raise RuntimeError(
            f"Could not find any SAR sources matching the provided parameters ({sar_input_params})")
    return sar_sources


def thread_download(
        downloaded_queue: multiprocessing.Queue,
        cleanup_queue: multiprocessing.Queue,
        sar_source: sar_data.SARSource,
        wind_source: Union[wind_data.ERA5Source, wind_data.Sentinel1OCNSource],
        input_dir: Path,
        wind_folder: Path,
        cleanup: bool = True):
    """Downloads SAR and wind. Meant to be run in a thread.
    """
    try:
        sar_source.download(input_dir, unzip=True)
    except requests.RequestException:
        logger.error("Error while downloading %s", sar_source.identifier, exc_info=True)
        if cleanup:
            cleanup_queue.put((sar_source, wind_source))
        return None

    try:
        wind_source.download(wind_folder)
    except RuntimeError:
        logger.warning("No wind data found matching %s", sar_source.identifier)
        if cleanup:
            cleanup_queue.put((sar_source, wind_source))
        return None

    downloaded_queue.put((sar_source, wind_source))


def multiprocess_preprocess(
        downloaded_queue: multiprocessing.Queue,
        preprocessed_queue: multiprocessing.Queue,
        cleanup_queue: multiprocessing.Queue,
        denoised_dir: Path,
        polarization: str,
        wind_folder: Path,
        pixel_size: int = 500,
        cleanup: bool = True):
    """Preprocess SAR (denoise, resize) and wind data (interpolate),
    then generate the full dataset.
    Meant to be run in a separate process
    """
    logger.debug("Starting preprocessing process")
    while True:
        next_item = downloaded_queue.get()
        if next_item is Done:
            logger.debug("Stopping preprocessing process")
            break
        sar_source, wind_source = next_item
        try:
            sar_source.preprocess(denoised_dir, polarization=polarization, pixel_size=pixel_size)
            wind_source.interpolate_on_sar_grid(wind_folder)
            preprocessed_queue.put((sar_source, wind_source))
        except Exception:
            logger.error("Error during preprocessing of %s", sar_source.identifier, exc_info=True)
            if cleanup:
                cleanup_queue.put((sar_source, wind_source))


def multiprocess_make_dataset(
        preprocessed_queue: multiprocessing.Queue,
        processed_queue: multiprocessing.Queue,
        cleanup_queue: multiprocessing.Queue,
        output_dir: Path,
        gmf: str,
        iterations: int,
        product_version: str = '1.0',
        file_version: str = '1.0',
        cleanup: bool = True,
        plot: bool = False):
    """Create the full dataset. Meant to be run in a separate process
    """
    logger.debug("Starting dataset making process")
    while True:
        next_item = preprocessed_queue.get()
        if next_item is Done:
            logger.debug("Stopping dataset making process")
            break
        sar_source, wind_source = next_item
        try:
            full_dataset_path = make_full_dataset(
                sar_source=sar_source, wind_source=wind_source,
                out_dir=output_dir,
                full_ds_file_name=make_product_file_name(sar_source, product_version, file_version),
                gmf=gmf, iterations=iterations)
            if plot:
                processed_queue.put((sar_source, wind_source, full_dataset_path))
            elif cleanup:
                cleanup_queue.put((sar_source, wind_source))
        except Exception:
            logger.error("Error during dataset creation for %s",
                         sar_source.identifier, exc_info=True)
            if cleanup:
                cleanup_queue.put((sar_source, wind_source))


def multiprocess_plot_dataset(
        processed_queue: multiprocessing.Queue,
        plot_dir: Path):
    """Plot the full dataset. Meant to be run in a separate process
    """
    logger.debug("Starting plotting process")
    while True:
        next_item = processed_queue.get()
        if next_item is Done:
            logger.debug("Stopping plotting process")
            break
        sar_source, wind_source, full_dataset_path = next_item
        try:
            plot_full_dataset(sar_source, wind_source, full_dataset_path, plot_dir)
        except Exception:
            logger.error("Error during plotting of %s", sar_source.identifier, exc_info=True)


def multiprocess_cleanup(cleanup_queue: multiprocessing.Queue):
    """Cleanup input and intermediary files after processing"""
    logger.debug("Starting cleanup process")
    while True:
        next_item = cleanup_queue.get()
        if next_item is Done:
            logger.debug("Stopping plotting process")
            break
        sar_source, wind_source = next_item
        logger.info("Deleting temporary files for %s", sar_source.identifier)
        to_delete: tuple[Path] = (
            sar_source.data_path,
            sar_source.data_path.parent / f"{sar_source.data_path.stem}.zip",
            sar_source.preprocessed_path,
            wind_source.data_path,
            wind_source.interpolated_path,
        )
        for p in to_delete:
            if p:
                logger.debug("Deleting %s", p)
                try:
                    p.unlink()
                except IsADirectoryError:
                    shutil.rmtree(p)
                except FileNotFoundError as e:
                    logger.debug("File does not exist, can't remove: %s", e.filename)
                except Exception:
                    logger.error("Error while deleting %s", p, exc_info=True)


def stop_processes(processes: dict, timeout: int = 1800):
    """Stop processes listenting to a queue"""
    # for each set of workers, send a message to stop and wait for them
    # to be finished before stopping the next workers
    for process_config in processes.values():
        # send stop messages in the input queue
        for _ in range(process_config['workers']):
            process_config['input_queue'].put(Done)
        # wait for the processes to stop
        for p in process_config['processes']:
            p.join(timeout)


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
        cleanup: bool = True,
        workdir: Union[str, Path] = Path('.'),
        input_dir: Union[str, Path] = None,
        output_dir: Union[str, Path] = None,
        denoised_dir: Union[str, Path] = None,
        wind_folder: Union[str, Path] = None,
        plot_dir: Union[str, Path] = None,
        max_download_workers: int = 10,
        max_preprocess_workers: int = 5,
        max_process_workers: int = 10,
        max_plot_workers: int = 5,
        pixel_size: int = 500,
        product_version: str = '1.0',
        file_version: str = '1.0'):
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

    downloaded_queue = multiprocessing.Queue()
    preprocessed_queue = multiprocessing.Queue()
    processed_queue = multiprocessing.Queue()
    cleanup_queue = multiprocessing.Queue()

    processes = {
        'preprocess': {
            'input_queue': downloaded_queue,
            'workers': max_preprocess_workers,
            'function': multiprocess_preprocess,
            'kwargs': {
                "downloaded_queue": downloaded_queue,
                "preprocessed_queue": preprocessed_queue,
                "cleanup_queue": cleanup_queue,
                "denoised_dir": denoised_dir,
                "polarization": 'VV',
                "wind_folder": wind_folder,
                "pixel_size": pixel_size,
                "cleanup": True,
            },
            'processes': [],
        },
        'process': {
            'input_queue': preprocessed_queue,
            'workers': max_process_workers,
            'function': multiprocess_make_dataset,
            'kwargs': {
                "preprocessed_queue": preprocessed_queue,
                "processed_queue": processed_queue,
                "cleanup_queue": cleanup_queue,
                "output_dir": output_dir,
                "gmf": gmf,
                "iterations": iterations,
                "product_version": product_version,
                "file_version": file_version,
                "cleanup": True,
                "plot": False,
            },
            'processes': [],
        },
    }

    if plot:
        processes['plot'] = {
            'input_queue': processed_queue,
            'workers': max_plot_workers,
            'function': multiprocess_plot_dataset,
            'kwargs': {
                "processed_queue": processed_queue,
                "plot_dir": plot_dir,
            },
            'processes': [],
        }
        # download coastlines
        cartopy.feature.download.__main__.download_features(['physical'])

    if cleanup:
        processes['cleanup'] = {
            'input_queue': cleanup_queue,
            'workers': 1,
            'function': multiprocess_cleanup,
            'kwargs': {"cleanup_queue": cleanup_queue},
            'processes': [],
        }

    with ThreadPoolExecutor(max_workers=max_download_workers) as download_executor:
        download_futures = []

        try:
            # start worker processes
            for process_config in processes.values():
                for _ in range(process_config['workers']):
                    p = multiprocessing.Process(
                        target=process_config['function'],
                        kwargs=process_config['kwargs'])
                    process_config['processes'].append(p)
                    p.start()

            # download SAR and wind, starting the processing chain
            for sar_source in sar_sources:
                full_path = check_final_file(
                    out_dir=output_dir,
                    sar_source=sar_source,
                    product_version=product_version,
                    file_version=file_version,
                    remove_invalid=True)
                if full_path:
                    continue

                wind_source = wind_source_class(sar_source)
                download_futures.append(download_executor.submit(
                    thread_download,
                    downloaded_queue=downloaded_queue,
                    cleanup_queue=cleanup_queue,
                    sar_source=sar_source, wind_source=wind_source,
                    input_dir=input_dir, wind_folder=wind_folder,
                    cleanup=cleanup))

            for download_future in as_completed(download_futures):
                try:
                    download_future.result()
                except Exception:
                    logger.error("Error during download", exc_info=True)

            # stop worker processes
            stop_processes(processes)

        except KeyboardInterrupt:
            for download_future in download_futures:
                download_future.cancel()
            stop_processes(processes)

    logger.info('Done processing sar sources matching %s', sar_input_params)
