#!/usr/bin/env python
# coding=utf-8
#-------------------------------------------------------------------------------
# Name:		cdop.py
# Purpose:      CDOP gmf to estimate the Doppler shift at C-band 
#
# Author:       Alexis Mouche
# Modified:	Morten Wergeland Hansen
#
# Created:	03.07.2014
# Last modified:03.02.2015 10:58
# Copyright:    (c) Ifremer
# License:      GNU General Public License, v.3
#               http://www.gnu.org/licenses/gpl-3.0.html
#-------------------------------------------------------------------------------

''' 
    This code contains the empirical CDOP geophysical model function, and is
    based on the following paper:
    
    Mouche A.A., Collard F., Chapron B., Dagestad K.-F., Guitton G.,
    Johannessen J.A., Kerbaol V., Hansen M.W. (2012). On the use of Doppler
    shift for sea surface wind retrieval from SAR. IEEE Transactions on
    Geoscience and Remote Sensing, Vol. 50, No. 7, pp 2901-2909,
    DOI:10.1109/TGRS.2011.2174998.

    CDOP estimates the Doppler shift at C-band given a wind field, the radar
    incidence angle and the radar transmit and receive polarizations (VV or
    HH).

'''
import numpy as np


def cdop_func(x):
    """
    """
    return 1./(1.+np.exp(-x))


def cdop(u10, phi, inc, pol):
    '''
        Input
        -----------
        u10 : float, numpy.array
              wind speed in m/s
        phi : float, numpy.array
              wind direction relative to SAR look direction
        inc : float, numpy.array
              SAR incidence angle
        pol : string
              SAR polarization (VV or HH)

        Return
        -----------
        dop : numpy.array
              Estimated Doppler shift [Hz] positive for surface motion toward the radar
    '''
    # Check inputs
    sizes = np.array([np.size(inc), np.size(u10), np.size(phi)])
    size = sizes.max()
    if ((sizes != size) & (sizes != 1)).any():
        raise Exception('Inputs sizes do not agree.')
    if pol.upper() not in ['VV', 'HH']:
        raise Exception('Unknown polarisation : '+pol)
    # NN coefficients (W=weights and B=biases)
    # (coefficient names in mouche2012 are given)
    if pol.upper() == 'VV':
        # lambda[0:2,1]
        B1 = np.array([-0.343935744939, 0.108823529412, 0.15],
                      dtype='float32')
        # lambda[0:2,0]
        W1 = np.array([0.028213254683, 0.0411764705882, .00388888888889],
                      dtype='float32')
        # omega[i,0]
        B2 = np.array([14.5077150927, -11.4312028555, 1.28692747109,
                       -1.19498666071, 1.778908726, 11.8880215573,
                       1.70176062351, 24.7941267067, -8.18756617111,
                       1.32555779345, -9.06560116738],
                      dtype='float32')
        # omega[i,[3,2,1]]
        W2 = np.array([[19.7873046673, 22.2237414308, 1.27887019276],
                       [2.910815875, -3.63395681095, 16.4242081101],
                       [1.03269004609, 0.403986575614, 0.325018607578],
                       [3.17100261168, 4.47461213024, 0.969975702316],
                       [-3.80611082432, -6.91334859293, -0.0162650756459],
                       [4.09854466913, -1.64290475596, -13.4031862615],
                       [0.484338480824, -1.30503436654, -6.04613303002],
                       [-11.1000239122, 15.993470129, 23.2186869807],
                       [-0.577883159569, 0.801977535733, 6.13874672206],
                       [0.61008842868, -0.5009830671, -4.42736737765],
                       [-1.94654022702, 1.31351068862, 8.94943709074]],
                      dtype='float32')
        # gamma[0]
        B3 = np.array(4.07777876994, dtype='float32')
        # gamma[1:11]
        W3 = np.array([7.34881153553, 0.487879873912, -22.167664703,
                       7.01176085914, 3.57021820094, -7.05653415486,
                       -8.82147148713, 5.35079872715, 93.627037987,
                       13.9420969201, -34.4032326496],
                      dtype='float32')
        # beta
        B4 = np.array(-52.2644487109, dtype='float32')
        # alpha
        W4 = np.array(111.528184073, dtype='float32')
    elif pol.upper() == 'HH':
        # lambda[0:2,1]
        B1 = np.array([-0.342097701547, 0.118181818182, 0.15],
                      dtype='float32')
        # lambda[0:2,0]
        W1 = np.array([0.0281843837385, 0.0318181818182, 0.00388888888889],
                      dtype='float32')
        # omega[i,0]
        B2 = np.array([1.30653883096, -2.77086154074, 10.6792861882,
                       -4.0429666906, -0.172201666743, 20.4895916824,
                       28.2856865516, -3.60143441597, -3.53935574111,
                       -2.11695768022, -2.57805898849],
                      dtype='float32')
        # omega[i,[3,2,1]]
        W2 = np.array([[-2.61087309812, -0.973599180956, -9.07176856257],
                       [-0.246776181361, 0.586523978839, -0.594867645776],
                       [17.9261562541, 12.9439063319, 16.9815377306],
                       [0.595882115891, 6.20098098757, -9.20238868219],
                       [-0.993509213443, 0.301856868548, -4.12397246171],
                       [15.0224985357, 17.643307099, 8.57886720397],
                       [13.1833641617, 20.6983195925, -15.1439734434],
                       [0.656338134446, 5.79854593024, -9.9811757434],
                       [0.122736690257, -5.67640781126, 11.9861607453],
                       [0.691577162612, 5.95289490539, -16.0530462],
                       [1.2664066483, 0.151056851685, 7.93435940581]],
                      dtype='float32')
        # gamma[0]
        B3 = np.array(2.68352095337, dtype='float32')
        # gamma[1:11]
        W3 = np.array([-8.21498722494, -94.9645431048, -17.7727420108,
                       -63.3536337981, 39.2450482271, -6.15275352542,
                       16.5337543167, 90.1967379935, -1.11346786284,
                       -17.57689699, 8.20219395141],
                      dtype='float32')
        # beta
        B4 = np.array(-66.9554922921, dtype='float32')
        # alpha
        W4 = np.array(136.216953823, dtype='float32')
    # Make inputs as a single matrix (and clip phi in [0,180])
    inputs = np.zeros((3, size), dtype='float32')
    for ivar, var in enumerate((inc, u10, phi)):
        if sizes[ivar] == 1:
            inputs[ivar, :] = np.repeat(var, size)
        else:
            inputs[ivar, :] = np.ravel(var)
        if ivar == 2:
            inputs[ivar, :] = np.abs(((inputs[ivar, :]+180) % 360)-180)
        inputs[ivar, :] *= W1[ivar]
        inputs[ivar, :] += B1[ivar]
    # Compute CDOP
    B2 = np.tile(B2.reshape((11, 1)), (1, size))
    dop = W4*cdop_func(np.dot(W3, cdop_func(np.dot(W2, inputs) + B2)) + B3) + B4
    # Reshape output
    # (using the shape of input which have the maximum ndim)
    ndims = np.array([np.ndim(inc), np.ndim(u10), np.ndim(phi)])
    tmp = np.where(sizes == size)[0]
    ivar = tmp[ndims[tmp].argmax()]
    shp = np.shape((inc, u10, phi)[ivar])
    dop = dop.reshape(shp)
    return dop

