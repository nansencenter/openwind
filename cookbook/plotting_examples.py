import numpy as np
import timeit
import matplotlib.pyplot as plt

from nansat.domain import Domain
from nansat.nsr import NSR
from openwind.sar_wind import SARWind
from nansat.nansatmap import Nansatmap
import openwind_integration_tests.openwind_test_archive as ota

test_data = ota.OpenWindTestData()

# Some simple time tests - should be cleaned up...
setup_small = '''
from openwind.sar_wind import SARWind
import openwind_integration_tests.openwind_test_archive as ota
test_data = ota.OpenWindTestData()
test_data.get_sentinel1a(fsize='small')
SARWind(test_data.sentinel1a['small'])
'''
setup_medium = '''
from openwind.sar_wind import SARWind
import openwind_integration_tests.openwind_test_archive as ota
test_data = ota.OpenWindTestData()
test_data.get_sentinel1a(fsize='medium')
SARWind(test_data.sentinel1a['medium'])
'''

def measure_time_example(fsize='small'):
    print(timeit.timeit(setup_small, number=10)/10.)

def measure_time_example(fsize='medium'):
    print(timeit.timeit(setup_medium, number=10)/10.)

###

def plot_s1a_example(fsize='small'):
    test_data.get_sentinel1a(fsize=fsize)
    w = SARWind(test_data.sentinel1a[fsize])
    cc = w.get_corners()
    lonmin = np.int(np.floor(np.min(cc[0])*100))/100.
    lonmax = np.int(np.ceil(np.max(cc[0])*100))/100.
    latmin = np.int(np.floor(np.min(cc[1])*100))/100.
    latmax = np.int(np.ceil(np.max(cc[1])*100))/100.
    w.reproject( Domain(NSR().wkt, ext='-lle %s %s %s %s -ts %s %s' %(lonmin,
        latmin, lonmax, latmax, (lonmax-lonmin)*110., (latmax-latmin)*110.) ) )
    u = w['U']
    v = w['V']
    nmap = Nansatmap(w, resolution='h')
    nmap.pcolormesh(np.hypot(u,v), vmax=18)
    nmap.add_colorbar(fontsize=8)
    nmap.quiver(u, v, step=20)#, scale=1, width=0.001)
    nmap.draw_continents()
    nmap.drawgrid()
    #nmap.drawmeridians(np.arange(lonmin, lonmax, 5), labels=[False, False,
    #    True, False])
    #nmap.drawparallels(np.arange(latmin, latmax, 3), labels=[True, False,
    #    False, False])

    # set size of the figure (inches)
    #nmap.fig.set_figheight(20)
    #nmap.fig.set_figwidth(15)

    # save figure to a PNG file
    nmap.draw_continents()
    #plt.suptitle(
    #    'High resolution\nwind speed and direction\nfrom Sentinel-1A and ' \
    #        'NCEP\n%s' %w.time_coverage_start.isoformat(),
    #    fontsize=8
    #)
    nmap.fig.savefig('s1a_wind_%s.png'%fsize, dpi=150, bbox_inches='tight')
    #nmap.drawgrid()
