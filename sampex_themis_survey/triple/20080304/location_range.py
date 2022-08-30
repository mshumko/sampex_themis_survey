"""
Plot the locations of the THEMIS ASIs, SAMPEX, and THEMIS probe footprints during 
the 2008-03-04 triple conjunction.
"""
from datetime import datetime

import asilib
import sampex
import matplotlib.pyplot as plt

from sampex_themis_survey.themis.footprint import THEMIS_footprint

themis_asis = ['FSMI']#, 'YKNF']
themis_probes = ['D', 'E', 'C']  # Also try 'B', 'C'
alt = 110  # km
# time_range = [datetime(2008, 3, 4, 4, 0, 0), datetime(2008, 3, 4, 6, 0, 0)]
time_range = [datetime(2008, 3, 4, 5, 49, 39), datetime(2008, 3, 4, 5, 55, 0)]

themis_footprints = {}

ax = asilib.make_map(land_color='grey', lat_bounds=(52, 70), lon_bounds=(-128, -90))

for themis_asi in themis_asis:
    t, _, _, _, _ = asilib.plot_map('THEMIS', themis_asi, time_range[0], alt, 
            ax=ax, asi_label=True, color_bounds=None, pcolormesh_kwargs={'rasterized':True})

for themis_probe in themis_probes:
    f = THEMIS_footprint(themis_probe, time_range)
    f.map(alt=alt, model='T89', mag_input={'Kp':1.667})
    themis_footprints[themis_probe] = {'time':f.time, 'footprint':f.footprint}
    ax.plot(f.footprint[:, 1], f.footprint[:, 0], 
        label=f'THEMIS-{themis_probe.upper()}',
        zorder=5)

plt.legend()
plt.show()