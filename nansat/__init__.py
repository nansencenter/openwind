from __future__ import absolute_import

from pkgutil import extend_path
__path__ = extend_path(__path__, __name__)

# repetition of code in nansat/nansat/__init__.py
from nansat.nsr import NSR
from nansat.domain import Domain
from nansat.nansat import Nansat

__all__ = ['NSR', 'Domain', 'Nansat']

