class TooHighResolutionError(Exception):
    def __init__(self, file='', msg=''):
        self.file = file
        self.msg = msg

class PolarizationError(Exception):
    def __init__(self, file='', msg=''):
        self.file = file
        self.msg = msg

