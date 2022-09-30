from datetime import datetime, timedelta


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
