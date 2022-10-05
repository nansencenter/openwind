from datetime import datetime, timedelta
from argparse import ArgumentParser


def measure_time(func):
    """Measure function execution time"""
    def wrapper(*args):
        start_time = datetime.now()
        result = func(*args)
        print(f'{datetime.now() - start_time} complete!')
        return result
    return wrapper


def check_inputs(func):
    def wrapper(**kwargs):
        # Check download arguments
        if 'download' in kwargs.keys():
            if kwargs['download'] and kwargs['dst'] is None:
                raise ValueError(f'Data download is set to True but the dst={kwargs["dst"]}')
        result = func(**kwargs)
        return result
    return wrapper


def round_time(timestamp: datetime) -> datetime:
    """
    Round datetime to the nearest hour

    :param timestamp: input datetime
    :returns timestamp: datetime rounded to the nearest hour
    """
    # If minutes are more then 30 then add an hour
    if timestamp.minute >= 30:
        timestamp = timestamp + timedelta(hours=1)
    # Set minutes seconds and microseconds to 0
    timestamp = timestamp.replace(microsecond=0, second=0, minute=0)
    return timestamp


def create_argparser() -> ArgumentParser:
    """
    Create a parser for inpus for the SAR wind processing

    :returns parser: command line argument parser
    """
    parser = ArgumentParser(
        description='Derive near-surface wind product from SAR backscatter and model wind' \
                    'product using openwind. The openwind relies on the Nansat for data ' \
                    'reading and reprojection. Hence it crucial that the profided files can ' \
                    'be read from using Nansat.')
    parser.add_argument('--sar_source', metavar='/path/to/sar/data', type=str, nargs='*', 
                        help='Path to SAR product. Can be a number of paths in which case each' \
                             'product will be processed separately')
    parser.add_argument('--wind_source', metavar='/path/to/wind/data', type=str, nargs='*', 
                        help='Path to the wind field at the local disc. If not provided ' \
                             'then ERA5 wind field will be automatically collocated and ' \
                             'used for retrieval. NOTE: Number of wind sources must be equival' \
                             'to the number of SAR files provided in <sar_source>')
    parser.add_argument('--pixel_size', metavar='N', type=float, nargs='?', default=1000,
                        help='Pixel size of the final product in meters. Note that the higher' \
                             'resolution significantly increases processing coasts (time/memory)' \
                             'If not provided, then default of 1000 m is used')
    parser.add_argument('--export_dst', metavar='/path/to/product/dir', type=str, nargs='?', default='/src',
                        help='Indicate storage directory for the product. If not provided then' \
                             'default of /src is used')
    return parser
