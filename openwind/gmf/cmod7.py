import math
import sys
import logging
from dataclasses import dataclass
from pathlib import Path

import numpy as np
from scipy.interpolate import RegularGridInterpolator


logger = logging.getLogger(__name__)


@dataclass
class LUTDimension():
    size: int
    min: float
    max: float
    step: float


def get_lut():
    """"""
    data_path = Path(__file__).parent / 'cmod7_data' / f"gmf_cmod7_vv.dat_{sys.byteorder}_endian"
    # Dimensions of GMF table
    m = 250  # wind speed min/max = 0.2-50 (step 0.2) [m/s] --> 250 pts
    n = 73   # dir min/max = 0-180 (step 2.5) [deg]   -->  73 pts
    p = 51   # inc min/max = 16-66 (step 1) [deg]     -->  51 pts
    gmf_table = np.fromfile(data_path, dtype=np.float32)
    # Remove head and tail
    gmf_table = gmf_table[1:-1]
    # To access the table as a three-dimensional Fortran-ordered m x n x p matrix,
    # reshape it
    gmf_table = gmf_table.reshape((m, n, p), order="F")
    return (gmf_table,
            LUTDimension(m, 0.2, 50., 0.2),
            LUTDimension(n, 0., 180., 2.5),
            LUTDimension(p, 16., 66., 1.))


def make_interpolator():
    """"""
    lut, speed_dim, dir_dim, inc_dim = get_lut()
    return RegularGridInterpolator(
        (np.arange(speed_dim.min, speed_dim.max + speed_dim.step, speed_dim.step),
         np.arange(dir_dim.min, dir_dim.max + dir_dim.step, dir_dim.step),
         np.arange(inc_dim.min, inc_dim.max + inc_dim.step, inc_dim.step)),
        lut,
        bounds_error=False,
        fill_value=None)


def cmod7_forward(speed, wind_dir, incidence, interpolator=None):
    """"""
    if interpolator is None:
        interpolator = make_interpolator()
    return interpolator((speed, wind_dir, incidence))


def cmod7_inverse(sigma0_obs, wind_dir, incidence, iterations=10):
    """"""
    interp = make_interpolator()

    # First guess wind speed
    V = np.array([10.]) * np.ones(sigma0_obs.shape)
    step=5

    # Iterating until error is smaller than threshold
    for iterno in range(1, iterations):
        # sigma0_calc = cmod7_forward(V, wind_dir, incidence, interpolator=interp)
        sigma0_calc = interpolate_3d(V, wind_dir, incidence)
        ind = sigma0_calc - sigma0_obs > 0
        V = V + step
        V[ind] = V[ind] - 2*step
        step = step/2

    return V


def interpolate_3d(speed, wind_dir, incidence):
    """Python implementation of the interpolation function used by KNMI Scatterometer Team.
    Defaults to boundary values for out of bound interpolation
    """
    lut, speed_dim, dir_dim, inc_dim = get_lut()

    p = (speed - speed_dim.min) / speed_dim.step
    q = (wind_dir - dir_dim.min) / dir_dim.step
    r = (incidence - inc_dim.min) / inc_dim.step

    i = np.floor(p).astype(int)
    j = np.floor(q).astype(int)
    k = np.floor(r).astype(int)

    for index, coeff, dim in ((i, p, speed_dim), (j, q, dir_dim), (k, r, inc_dim)):
        index_negative_mask = index < 0
        index[index_negative_mask] = 0
        coeff[index_negative_mask] = 0.0
        index[index >= dim.size - 1] = dim.size - 2
        coeff[~index_negative_mask] -= index[~index_negative_mask]

    sigma111 = lut[i, j, k]
    sigma211 = lut[i+1, j, k]
    sigma121 = lut[i, j+1, k]
    sigma112 = lut[i, j, k+1]
    sigma221 = lut[i+1 , j+1, k]
    sigma212 = lut[i+1, j, k+1]
    sigma122 = lut[i, j+1, k+1]
    sigma222 = lut[i+1, j+1, k+1]

    sigma0 = (sigma111 +
              p*(sigma211-sigma111) +
              q*(sigma121-sigma111) +
              r*(sigma112-sigma111) +
              p*q*(sigma221+sigma111-sigma121-sigma211) +
              p*r*(sigma212+sigma111-sigma112-sigma211) +
              q*r*(sigma122+sigma111-sigma112-sigma121) +
              p*q*r*(sigma222+sigma112+sigma121+sigma211 -
                     sigma122-sigma212-sigma221-sigma111))

    return sigma0
