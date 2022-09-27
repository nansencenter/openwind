from datetime import datetime


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