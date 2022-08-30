"""
Plot the locations of the THEMIS ASIs, SAMPEX, and THEMIS probe footprints during 
the 2008-03-04 triple conjunction.
"""
from datetime import datetime

import asilib
import sampex
import matplotlib.pyplot as plt

from sampex_themis_survey.themis.footprint import THEMIS_footprint

themis_asis = ['FSMI', 'YKNF']
themis_probes = ['D', 'E']  # And 'B', 'C'
time_range = [datetime(2008, 3, 4, 4, 0, 0), datetime(2008, 3, 4, 6, 0, 0)]

ax = asilib.make_map(land_color='w')

for themis_probe in themis_probes:
    f = THEMIS_footprint(themis_probe, time_range)
    f.map()
    ax.plot(f.footprint[:, 1], f.footprint[:, 0], label=f'THEMIS-{themis_probe.upper()}')

plt.show()