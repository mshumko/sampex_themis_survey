"""
Plot the locations of the THEMIS ASIs, SAMPEX, and THEMIS probe footprints during 
the 2008-03-04 triple conjunction.
"""
from datetime import datetime

import sampex
from sampex_themis_survey.themis.footprint import THEMIS_footprint

themis_asis = ['FSMI', 'YKNF']
themis_probes = ['D', 'E']  # And 'B', 'C'
time_range = [datetime(2008, 3, 4, 4, 0, 0), datetime(2008, 3, 4, 6, 0, 0)]

for themis_probe in themis_probes:
    f = THEMIS_footprint(themis_probe, time_range)
    f.map()
    pass