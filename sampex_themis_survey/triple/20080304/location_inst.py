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
import matplotlib.colors
import matplotlib.dates

from sampex_themis_survey.themis.footprint import THEMIS_footprint
from sampex_themis_survey.themis.sst import SST
from sampex_themis_survey.themis.fbk import FBK

themis_asi_locations = ['FSMI']#, 'YKNF']
themis_probes = ['C'] #['D', 'E', 'C']  # Also try 'B', 'C'
alt = 110  # km
times = [
    datetime(2008, 3, 4, 5, 30, 0),
    datetime(2008, 3, 4, 5, 40, 0),
    datetime(2008, 3, 4, 5, 50, 0),
    datetime(2008, 3, 4, 6, 0, 0)
    ]
plot_range = (datetime(2008, 3, 4, 5, 0, 0), datetime(2008, 3, 4, 6, 30, 0))
lat_bounds=(55, 66)
lon_bounds=(-121, -101)

n = len(times)
fig = plt.figure(figsize=(12, 9))
spec = gridspec.GridSpec(
    nrows=len(themis_asi_locations)+4, 
    ncols=n, figure=fig, height_ratios=(2, 1, 1, 1, 1)
    )

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
    
    ax[i].text(0, 0, f'({subplot_letter}) {image_time.strftime("%H:%M:%S")}', 
        transform=ax[i].transAxes, va='bottom', color='white', fontsize=15)

themis_footprints = {}
colors = ['r', 'c', 'b', 'k']
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
    bx.tick_params(axis="x", labelbottom=False)
    asilib.plot_keogram('THEMIS', themis_asi_location, plot_range, 
        ax=bx, map_alt=alt, aacgm=True)
    bx.set_xlim(times[0], times[-1])
    bx.set_title('')
    bx.text(0, 1, f'({string.ascii_uppercase[n+i-1]}) THEMIS ASI keogram', 
        transform=bx.transAxes, va='top', color='white', fontsize=15)
    bx.set_ylabel('Magnetic lat [geg]')
    bx.set_ylim(66, 72)

cx = [fig.add_subplot(spec[i, :], sharex=bx) 
    for i in range(2, len(themis_asi_locations)+4)]
for cx_i in cx[:-1]:
    cx_i.tick_params(axis="x", labelbottom=False) 

# Electrons
cx_c = cx[0].inset_axes([1.01, 0, 0.02, 1], transform=cx[0].transAxes)
s = SST(themis_probes[0], plot_range[0], species='e')
s.load()
_, p = s.spectrum(ax=cx[0], 
    pcolormesh_kwargs={'norm':matplotlib.colors.LogNorm(vmin=1E2, vmax=5E5)}
    )
plt.colorbar(p, ax=cx, cax=cx_c, label='Electron flux')
cx[0].set_xlim(plot_range)
cx[0].set_ylabel('Energy [keV]')
cx[0].set_yscale('log')

# Protons
cx_c = cx[1].inset_axes([1.01, 0, 0.02, 1], transform=cx[1].transAxes)
s = SST(themis_probes[0], plot_range[0], species='i')
s.load()
_, p = s.spectrum(ax=cx[1], 
    pcolormesh_kwargs={'norm':matplotlib.colors.LogNorm(vmin=0.1, vmax=1E6)}
    )
plt.colorbar(p, ax=cx[1], cax=cx_c, label='Proton flux')
cx[1].set_xlim(plot_range)
cx[1].set_ylabel('Energy [keV]')
cx[1].set_yscale('log')

# Waves
s = FBK(themis_probes[0], plot_range[0])
s.load()
s.spectrum(ax=cx[2], 
    pcolormesh_kwargs={'norm':matplotlib.colors.LogNorm(vmin=1E-4, vmax=0.1)})
cx[2].set_xlabel('Time [HH:MM]')
cx[2].set_ylabel('Frequency [Hz]') 
cx[2].set_yscale('log')
cx[2].set_ylim(1, 1E3)
time_format = matplotlib.dates.DateFormatter('%H:%M')
cx[2].xaxis.set_major_formatter(time_format)

plot_labels = (
    'THEMIS-SST Electrons',
    'THEMIS-SST Ions',
    'THEMIS-FBK'
)
themis_label_y = (0, 0, 0.95)

for i, (cx_i, plot_label, y) in enumerate(zip(cx, plot_labels, themis_label_y), start=n+1):
    cx_i.text(0, y, f'({string.ascii_uppercase[i]}) {plot_label}', 
        transform=cx_i.transAxes, va='bottom', color='white', fontsize=15)

# cx = fig.add_subplot(spec[-3, :], sharex=bx)
# # cx.get_shared_x_axes().join(bx, cx)
# dx = cx.inset_axes([1.01, 0, 0.02, 1], transform=cx.transAxes)
# s = SST(themis_probes[0], plot_range[0], species='e')
# s.load()
# _, p = s.spectrum(ax=cx, 
#     pcolormesh_kwargs={'norm':matplotlib.colors.LogNorm(vmin=1E2, vmax=5E5)}
#     )
# plt.colorbar(p, ax=cx, cax=dx, label='Electron flux')
# cx.set_xlim(plot_range)
# cx.set_ylabel('Energy [keV]')
# cx.set_yscale('log')
# cx.text(0, 0.8, f'({string.ascii_uppercase[n+1]}) THEMIS-SST Electrons', 
#         transform=cx.transAxes, va='top', color='white', fontsize=15)
# cx.tick_params(axis="x", labelbottom=False)      

# ex = fig.add_subplot(spec[-3, :], sharex=cx)
# # cx.get_shared_x_axes().join(bx, cx)
# fx = ex.inset_axes([1.01, 0, 0.02, 1], transform=ex.transAxes)
# s = SST(themis_probes[0], plot_range[0], species='i')
# s.load()
# _, p = s.spectrum(ax=fx, 
#     pcolormesh_kwargs={'norm':matplotlib.colors.LogNorm(vmin=0.1, vmax=1E6)}
#     )
# plt.colorbar(p, ax=fx, cax=ex, label='Proton flux')
# fx.set_xlim(plot_range)
# fx.set_ylabel('Energy [keV]')
# fx.set_yscale('log')
# fx.text(0, 0.8, f'({string.ascii_uppercase[n+2]}) THEMIS-SST Protons', 
#         transform=cx.transAxes, va='top', color='white', fontsize=15)
# cx.tick_params(axis="x", labelbottom=False)    

# fx = fig.add_subplot(spec[-1, :], sharex=cx)
# s = FBK(themis_probes[0], plot_range[0])
# s.load()
# s.spectrum(ax=fx, 
#     pcolormesh_kwargs={'norm':matplotlib.colors.LogNorm(vmin=1E-4, vmax=0.1)})
# fx.text(0, 0.8, f'({string.ascii_uppercase[n+2]}) THEMIS-FBK', 
#         transform=fx.transAxes, va='top', color='white', fontsize=15)
# fx.set_xlabel('Time [HH:MM]')
# fx.set_ylabel('Frequency [Hz]') 
# fx.set_yscale('log')
# fx.set_ylim(1, 1E3)
# time_format = matplotlib.dates.DateFormatter('%H:%M')
# fx.xaxis.set_major_formatter(time_format)

plt.suptitle(f'Example Triple Conjunction | THEMIS-THEMIS ASI-SAMPEX', fontsize=15)

# plt.legend()
plt.subplots_adjust(wspace=0.02, hspace=0.07, left=0.055, right=0.92, top=0.943)
plt.show()