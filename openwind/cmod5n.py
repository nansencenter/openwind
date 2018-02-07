from numpy import cos, exp, tanh, ones, array
import warnings
# Ignore overflow errors for wind calculations over land
warnings.simplefilter("ignore", RuntimeWarning)

# CMOD5.N model coefficients
# https://www.ecmwf.int/en/elibrary/9873-cmod5n-c-band-geophysical-model-function-equivalent-neutral-wind
# NB: 0 added as first element below, to avoid switching from 1-indexing to 0-indexing
C = [0, -0.6878, -0.7957, 0.3380, -0.1728, 0.0000, 0.0040, 0.1103, 0.0159, 6.7329, 2.7713,
     -2.2885, 0.4971, -0.7250, 0.0450, 0.0066, 0.3222, 0.0120, 22.7000, 2.0813, 3.0000,
     8.3659, -3.3428, 1.3236, 6.2437, 2.3893, 0.3249, 4.1590, 1.6930]


def cmod5n_forward(v, phi, theta):
    """cmod5n_forward(v, phi, theta)
    Parameters
    ----------
        v: float, numpy.array
            2D array with wind velocities in [m/s] (always >= 0)
        phi: float, numpy.array
            2D array with angles between azimuth and wind direction in [deg] (= D - AZM)
        theta: float, numpy.array
            2D array with incidence angles in [deg]

    Returns
    -------
        cmod_n: float, numpy.array
            2D array with normalized backscatter (linear)

    The CMOD5 Model Formulation and Coefficients were published in:
    H. Hersbach, A. Stoffelen, and S. de Haan. 2007. An improved C-band scatterometer ocean
    geophysical model function: CMOD5. Journal of Geophysical Research, Vol. 112, C03006,
    doi:10.1029/2006JC003743

    A. STOFFELEN              MAY  1991 ECMWF  CMOD4
    A. STOFFELEN, S. DE HAAN  DEC  2001 KNMI   CMOD5 PROTOTYPE
    H. HERSBACH               JUNE 2002 ECMWF  COMPLETE REVISION
    J. de Kloe                JULI 2003 KNMI,  rewritten in fortan90
    A. Verhoef                JAN  2008 KNMI,  CMOD5 for neutral winds
    K.F.Dagestad              OCT 2011 NERSC,  Vectorized Python version
    """

    DTOR = 57.29577951
    THETM = 40.
    THETHR = 25.
    ZPOW = 1.6

    Y0 = C[19]
    PN = C[20]

    A = C[19] - (C[19] - 1) / C[20]
    B = 1. / (C[20] * (C[19] - 1.) ** (3-1))

    # ANGLES
    FI = phi/DTOR
    CSFI = cos(FI)
    CS2FI = 2.0 * CSFI * CSFI - 1.0

    X = (theta - THETM) / THETHR

    # B0: FUNCTION OF WIND SPEED AND INCIDENCE ANGLE
    A0 = C[1] + C[2] * X + C[3] * (X ** 2) + C[4] * (X ** 3)
    A1 = C[5] + C[6] * X
    A2 = C[7] + C[8] * X

    GAM = C[9] + C[10] * X + C[11] * (X ** 2)
    S0 = C[12] + C[13] * X
    
    # V is missing! Using V=v as substitute, this is apparently correct
    V = v
    S = A2 * V
    S_vec = S.copy() 
    SlS0 = [S_vec < S0]
    S_vec[SlS0] = S0[SlS0]
    A3 = 1. / (1. + exp(-S_vec))
    SlS0 = (S < S0)
    A3[SlS0] = A3[SlS0] * (S[SlS0] / S0[SlS0]) ** (S0[SlS0] * (1. - A3[SlS0]))
    B0 = (A3 ** GAM) * 10. ** (A0 + A1 * V)
        
    # B1: FUNCTION OF WIND SPEED AND INCIDENCE ANGLE
    B1 = C[15] * V * (0.5 + X - tanh(4. * (X + C[16] + C[17] * V)))
    B1 = C[14] * (1. + X) - B1
    B1 = B1 / (exp(0.34 * (V - C[18])) + 1.)

    # B2: FUNCTION OF WIND SPEED AND INCIDENCE ANGLE
    V0 = C[21] + C[22] * X + C[23] * (X ** 2)
    D1 = C[24] + C[25] * X + C[26] * (X ** 2)
    D2 = C[27] + C[28] * X

    V2 = (V / V0 + 1.)
    V2ltY0 = V2 < Y0
    V2[V2ltY0] = A + B * (V2[V2ltY0] - 1.) ** PN
    B2 = (-D1 + D2 * V2) * exp(-V2)

    # CMOD5_N: COMBINE THE THREE FOURIER TERMS
    cmod5_n = B0 * (1.0 + B1 * CSFI + B2 * CS2FI) ** ZPOW

    return cmod5_n
    

def cmod5n_inverse(sigma0_obs, phi, incidence, iterations=10):
    """The function iterates the forward CMOD5N <cmod5n_forward> function until agreement with input
     (observed) sigma0 values

    Parameters:
        sigma0_obs: float, numpy.array
             2D array with Normalized Radar Cross Section (NRCS) [linear units]
        phi: float, numpy.array
            2D array with angles between azimuth and wind direction in [deg] (= D - AZM)
        incidence: float, numpy.array
            incidence angles in [deg]
        iterations: int
            number of iterations to run

    Returns:
        v: float, numpy.array
            2D array with wind speeds at 10 m, neutral stratification
    """

    # First guess wind speed
    v = array([10.]) * ones(sigma0_obs.shape)
    step = 10.
    
    # Iterating until error is smaller than threshold
    for iterno in range(1, iterations):
        #print iterno
        sigma0_calc = cmod5n_forward(v, phi, incidence)
        ind = sigma0_calc - sigma0_obs > 0
        v = v + step
        v[ind] = v[ind] - 2 * step
        step = step / 2

    return v
