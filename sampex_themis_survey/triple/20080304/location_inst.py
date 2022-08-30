"""
Plot the instantaneous locations of the THEMIS ASIs, and THEMIS probe footprints 
during  the 2008-03-04 triple conjunction.
"""
from datetime import datetime
import string

import numpy as np
import asilib
import sampex
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec

from sampex_themis_survey.themis.footprint import THEMIS_footprint

themis_asi_locations = ['FSMI']#, 'YKNF']
themis_probes = ['D', 'E', 'C']  # Also try 'B', 'C'
alt = 110  # km
times = [
    datetime(2008, 3, 4, 5, 30, 0),
    datetime(2008, 3, 4, 5, 40, 0),
    datetime(2008, 3, 4, 5, 50, 0),
    datetime(2008, 3, 4, 6, 0, 0)
    ]
lat_bounds=(52, 70)
lon_bounds=(-128, -90)

n = len(times)
fig = plt.figure(figsize=(12, 5))
spec = gridspec.GridSpec(nrows=len(themis_asi_locations)+1, ncols=n, figure=fig, height_ratios=(2, 1))

ax = n*[None]
nearest_asi_image_times = []
z = zip(times, string.ascii_uppercase[:n])

for i, (image_time, subplot_letter) in enumerate(z):
    ax[i] = fig.add_subplot(spec[0, i])
    ax[i].get_xaxis().set_visible(False)
    ax[i].get_yaxis().set_visible(False)
    asilib.make_map(lat_bounds=lat_bounds, lon_bounds=lon_bounds, ax=ax[i], land_color='grey')
    for themis_asi_location in themis_asi_locations:
        t, _, _, _, _ = asilib.plot_map('THEMIS', themis_asi_location, image_time, alt, 
            ax=ax[i], asi_label=False, color_bounds=None, pcolormesh_kwargs={'rasterized':True})
        if themis_asi_location == themis_asi_locations[0]:
            nearest_asi_image_times.append(t)
    
    ax[i].text(0, 1, f'({subplot_letter}) {image_time.strftime("%H:%M:%S")}', 
        transform=ax[i].transAxes, va='top', color='green', fontsize=15)

themis_footprints = {}
colors = ['c', 'g', 'b', 'k']
for themis_probe, color in zip(themis_probes, colors):
    f = THEMIS_footprint(themis_probe, (times[0], times[-1]))
    f.map(alt=alt, model='T89', mag_input={'Kp':1.667})
    themis_footprints[themis_probe] = {'time':f.time, 'footprint':f.footprint}
    
    for time, ax_i in zip(times, ax):
        dt = np.abs([(time-ti).total_seconds() for ti in f.time])
        idt = np.argmin(dt)
        assert (time-f.time[idt]).total_seconds() <= 60, (
            f'The time stamp difference, {(time-f.time[idt]).total_seconds()} '
            f's, is too large.')
        ax_i.scatter(f.footprint[idt, 1], f.footprint[idt, 0], 
            label=f'THEMIS-{themis_probe.upper()}',
            color=color, s=100, marker='.')

for i, themis_asi_location in enumerate(themis_asi_locations, start=1):
    bx = fig.add_subplot(spec[i, :])
    asilib.plot_keogram('THEMIS', themis_asi_location, (times[0], times[-1]), 
        ax=bx, map_alt=alt, aacgm=True)
    bx.set_xlim(times[0], times[-1])
    bx.set_title('')

plt.suptitle(f'Example Triple Conjunction | THEMIS probes | THEMIS ASI | SAMPEX', fontsize=15)

# plt.legend()
plt.show()