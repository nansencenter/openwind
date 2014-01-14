import warnings
# Ignore overflow errors for wind calculations over land
warnings.simplefilter("ignore", RuntimeWarning) 

def cmod5n_forward(v,phi,theta):
    '''!     ---------
    !     cmod5n_forward(v, phi, theta)
    !         inputs:
    !              v     in [m/s] wind velocity (always >= 0)
    !              phi   in [deg] angle between azimuth and wind direction
    !                    (= D - AZM)
    !              theta in [deg] incidence angle
    !         output:
    !              CMOD5_N NORMALIZED BACKSCATTER (LINEAR)
    !
    !        All inputs must be Numpy arrays of equal sizes
    !
    !     A. STOFFELEN              MAY  1991 ECMWF  CMOD4
    !     A. STOFFELEN, S. DE HAAN  DEC  2001 KNMI   CMOD5 PROTOTYPE
    !     H. HERSBACH               JUNE 2002 ECMWF  COMPLETE REVISION
    !     J. de Kloe                JULI 2003 KNMI,  rewritten in fortan90
    !     A. Verhoef                JAN  2008 KNMI,  CMOD5 for neutral winds
    !     K.F.Dagestad              OCT 2011 NERSC,  Vectorized Python version
    !---------------------------------------------------------------------
       '''
        
    from numpy import cos, exp, tanh, array
    
    DTOR   = 57.29577951
    THETM  = 40.
    THETHR = 25.
    ZPOW   = 1.6
    
    # NB: 0 added as first element below, to avoid switching from 1-indexing to 0-indexing
    C = [0, -0.6878, -0.7957,  0.3380, -0.1728, 0.0000,  0.0040, 0.1103, 0.0159, 
          6.7329,  2.7713, -2.2885, 0.4971, -0.7250, 0.0450, 
          0.0066,  0.3222,  0.0120, 22.7000, 2.0813,  3.0000, 8.3659,
          -3.3428,  1.3236,  6.2437,  2.3893, 0.3249,  4.1590, 1.6930]
    Y0 = C[19]
    PN = C[20]
    A  = C[19]-(C[19]-1)/C[20]

    B  = 1./(C[20]*(C[19]-1.)**(3-1))

#  !  ANGLES
    FI=phi/DTOR
    CSFI = cos(FI)
    CS2FI= 2.00 * CSFI * CSFI - 1.00

    X  = (theta - THETM) / THETHR
    XX = X*X

    #  ! B0: FUNCTION OF WIND SPEED AND INCIDENCE ANGLE
    A0 =C[ 1]+C[ 2]*X+C[ 3]*XX+C[ 4]*X*XX
    A1 =C[ 5]+C[ 6]*X
    A2 =C[ 7]+C[ 8]*X

    GAM=C[ 9]+C[10]*X+C[11]*XX
    S0 =C[12]+C[13]*X
    
    # V is missing! Using V=v as substitute, this is apparently correct
    V=v
    S = A2*V
    S_vec = S.copy() 
    SlS0 = [S_vec<S0]
    S_vec[SlS0]=S0[SlS0]
    A3=1./(1.+exp(-S_vec))
    SlS0 = (S<S0)
    A3[SlS0]=A3[SlS0]*(S[SlS0]/S0[SlS0])**( S0[SlS0]*(1.- A3[SlS0]))
    #A3=A3*(S/S0)**( S0*(1.- A3))
    B0=(A3**GAM)*10.**(A0+A1*V)
        
    #  !  B1: FUNCTION OF WIND SPEED AND INCIDENCE ANGLE
    B1 = C[15]*V*(0.5+X-tanh(4.*(X+C[16]+C[17]*V)))
    B1 = C[14]*(1.+X)- B1
    B1 = B1/(exp( 0.34*(V-C[18]) )+1.)

    #  !  B2: FUNCTION OF WIND SPEED AND INCIDENCE ANGLE
    V0 = C[21] + C[22]*X + C[23]*XX
    D1 = C[24] + C[25]*X + C[26]*XX
    D2 = C[27] + C[28]*X

    V2 = (V/V0+1.)
    V2ltY0 = V2<Y0
    V2[V2ltY0] = A+B*(V2[V2ltY0]-1.)**PN
    B2 = (-D1+D2*V2)*exp(-V2)

    #  !  CMOD5_N: COMBINE THE THREE FOURIER TERMS
    CMOD5_N = B0*(1.0+B1*CSFI+B2*CS2FI)**ZPOW
    return CMOD5_N
    

def cmod5n_inverse(sigma0_obs, phi, incidence, iterations=10):
    '''!     ---------
    !     cmod5n_inverse(sigma0_obs, phi, incidence, iterations)
    !         inputs:
    !              sigma0_obs     Normalized Radar Cross Section [linear units]
    !              phi   in [deg] angle between azimuth and wind direction
    !                    (= D - AZM)
    !              incidence in [deg] incidence angle
    !              iterations: number of iterations to run
    !         output:
    !              Wind speed, 10 m, neutral stratification 
    !
    !        All inputs must be Numpy arrays of equal sizes
    !
    !    This function iterates the forward CMOD5N function
    !    until agreement with input (observed) sigma0 values   
    !---------------------------------------------------------------------
       '''
    from numpy import ones, array
    
    # First guess wind speed
    V = array([10.])*ones(sigma0_obs.shape);
    step=10.
    
    # Iterating until error is smaller than threshold
    for iterno in range(1, iterations):
        #print iterno
        sigma0_calc = cmod5n_forward(V, phi, incidence)
        ind = sigma0_calc-sigma0_obs>0
        V = V + step
        V[ind] = V[ind] - 2*step 
        step = step/2

    #mdict={'s0obs':sigma0_obs,'s0calc':sigma0_calc}
    #from scipy.io import savemat
    #savemat('s0test',mdict)
        
    return V
