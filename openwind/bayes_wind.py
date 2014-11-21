#-------------------------------------------------------------------------------
# Name:		bayes_wind.py
# Purpose:      
#
# Author:       Morten Wergeland Hansen
# Modified:	Morten Wergeland Hansen
#
# Created:	21.11.2014
# Last modified:21.11.2014 17:04
# Copyright:    (c) NERSC
# License:      
#-------------------------------------------------------------------------------
import numpy as np
from openwind.sar_wind import SARWind

def cost_function(dom, obs, err):
    return np.square(( dom - obs )/err)


class BayesianWind(SARWind):
    wind_speed_domain = np.linspace(-20,20,401)
    model_wind_uncertainty = 1.5

    def __init__(self):
        [self.u_domain,self.v_domain] = np.meshgrid(wind_speed_domain,
                wind_speed_domain)

