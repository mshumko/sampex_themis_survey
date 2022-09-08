from datetime import datetime
import string

import matplotlib.pyplot as plt
import matplotlib.dates
from matplotlib.ticker import FuncFormatter
import matplotlib.patches
import matplotlib.gridspec as gridspec
import matplotlib.collections
import matplotlib.colors
import numpy as np
import pandas as pd
import asilib
import sampex

from sampex_themis_survey.footprint import SAMPEX_footprint

sampex_time_range = (datetime(2012, 3, 16, 9, 31, 0), datetime(2012, 3, 16, 9, 34, 15))
plot_time_range = (datetime(2012, 3, 16, 9, 31, 30), datetime(2012, 3, 16, 9, 33, 30))
sampex_x_labels = {'L':'L_Shell', 'MLT':'MLT', 'Geo Lat':'GEO_Lat', 'Geo Lon':'GEO_Long'}

color_footprint = False
image_times = [
            datetime(2012, 3, 16, 9, 31, 45),
            datetime(2012, 3, 16, 9, 32, 10),
            datetime(2012, 3, 16, 9, 32, 45)
            ]
n = len(image_times)
asi_array_code = 'THEMIS'
themis_location_code = 'KAPU'
map_alt = 110
lon_bounds = (-90, -71)
lat_bounds = (45, 61)
color_bounds = (3500, 4000)

fig = plt.figure(figsize=(10, 5))
spec = gridspec.GridSpec(nrows=2, ncols=n, figure=fig, height_ratios=(2, 1))

ax = n*[None]
nearest_asi_image_times = []
z = zip(ax, image_times, string.ascii_uppercase[:n])

for i, (ax_i, image_time, subplot_letter) in enumerate(z):
    ax[i] = fig.add_subplot(spec[0, i])
    ax[i].get_xaxis().set_visible(False)
    ax[i].get_yaxis().set_visible(False)
    asilib.make_map(lat_bounds=lat_bounds, lon_bounds=lon_bounds, ax=ax[i])
    t, _, _, _, _ = asilib.plot_map('THEMIS', themis_location_code, image_time, map_alt, 
        ax=ax[i], asi_label=False, color_bounds=None, pcolormesh_kwargs={'rasterized':True})
    asilib.plot_map('THEMIS', 'SNKQ', image_time, map_alt, 
        ax=ax[i], asi_label=False, color_bounds=None, pcolormesh_kwargs={'rasterized':True})
    
    nearest_asi_image_times.append(t)


# Now load, map, and interpolate the SAMPEX attitude data
asi_times, _ = asilib.load_image('THEMIS', themis_location_code, 
    time_range=(sampex_time_range[0]-pd.Timedelta(seconds=3), sampex_time_range[1]))
asi_times_df = pd.DataFrame(index=asi_times)
footprint = SAMPEX_footprint(sampex_time_range[0]).map(alt=map_alt)

# The merge will double the time cadence to 3 seconds and the new times will have
# an associated NaNs. df.interpolate will interpolate over the NaN values.
footprint = pd.merge_asof(
    asi_times_df, footprint, 
    left_index=True, right_index=True, 
    tolerance=pd.Timedelta(seconds=1), 
    direction='nearest'
    )
footprint = footprint.interpolate(method='time', limit_area='inside').dropna()

hilt = sampex.HILT(sampex_time_range[0]).load()
hilt = hilt.loc[sampex_time_range[0]:sampex_time_range[1], :]

def format_fn(tick_val, tick_pos):
    """
    The tick magic happens here. pyplot gives it a tick time, and this function 
    returns the closest label to that time. Read docs for FuncFormatter().
    """
    # Find the nearest time within 6 seconds (the cadence of the SAMPEX attitude files)
    tick_time = matplotlib.dates.num2date(tick_val).replace(tzinfo=None)
    i_min_time = np.argmin(np.abs(footprint.index - tick_time))
    if np.abs(footprint.index[i_min_time] - tick_time).total_seconds() > 3:
        # return ''
        raise ValueError(f'Nearest timestamp to tick_time is more than 6 seconds away')
    pd_index = footprint.index[i_min_time]
    # Cast np.array as strings so that it can insert the time string. 
    values = footprint.loc[pd_index, sampex_x_labels.values()].to_numpy().round(2).astype(str)
    values = np.insert(values, 0, tick_time.strftime('%H:%M:%S'))
    label = '\n'.join(values)
    return label

bx = fig.add_subplot(spec[-1, :])
bx.set_xlim(sampex_time_range)
# bx.set_yscale('log')
bx.set_ylim(9, 10**4)
bx.set_ylabel(f'>1 MeV electrons\n[counts/20 ms]')
bx.xaxis.set_major_formatter(FuncFormatter(format_fn))
bx.xaxis.set_major_locator(matplotlib.dates.SecondLocator(interval=10))

bx.set_xlabel('\n'.join(['Time'] + list(sampex_x_labels.keys())))
bx.xaxis.set_label_coords(-0.07, -0.09)

bx.text(0, 0.99, f'({string.ascii_uppercase[n]}) SAMPEX-HILT', va='top', 
    transform=bx.transAxes, weight='bold', fontsize=15)

if color_footprint:
    # Need to smooth so the microburst intervals are well defined.
    smooth_sec = 5
    smoothed_hilt = hilt.rolling(int(smooth_sec//20E-3), center=True).max()
    resampled_counts = pd.merge_asof(
        footprint, smoothed_hilt, 
        left_index=True, right_index=True, 
        tolerance=pd.Timedelta(seconds=20E-3), 
        direction='nearest'
        )
    # See example:
    # https://matplotlib.org/stable/gallery/lines_bars_and_markers/multicolored_line.html
    points = np.array(
        [resampled_counts.loc[:, 'GEO_Long'], resampled_counts.loc[:, 'GEO_Lat']]
        ).T.reshape(-1, 1, 2)
    segments = np.concatenate([points[:-1], points[1:]], axis=1)
    cmap = matplotlib.colors.ListedColormap(['r', 'g', 'b'])
    norm = matplotlib.colors.BoundaryNorm([1E2, 1e3], cmap.N)

# First row
for (ax_i, t, subplot_letter) in zip(ax, nearest_asi_image_times, string.ascii_uppercase[:n]):
    # Determine how to plot the entire footprint.
    if color_footprint:
        lc = matplotlib.collections.LineCollection(segments, cmap=cmap, norm=norm)
        lc.set_array(resampled_counts['counts'])
        lc.set_linewidth(2)
        line = ax_i.add_collection(lc)
    else:
        ax_i.plot(footprint.loc[:, 'GEO_Long'], footprint.loc[:, 'GEO_Lat'], 'r:')
    ax_i.scatter(
        footprint.loc[t, 'GEO_Long'], 
        footprint.loc[t, 'GEO_Lat'],
        c='red', s=150, marker='.',
        )

    ax_i.text(0, 1, f'({subplot_letter}) {t.strftime("%H:%M:%S")}', 
        transform=ax_i.transAxes, va='top', color='white', fontsize=15)

# Second row
bx.plot(hilt.index, hilt['counts'], c='r')
bx.xaxis.set_minor_locator(matplotlib.dates.SecondLocator())
bx.set_xlim(*plot_time_range)

for ax_i, image_time_numeric in zip(ax, matplotlib.dates.date2num(nearest_asi_image_times)):
    # Connecting lines between subplots
    line = matplotlib.patches.ConnectionPatch(
        xyA=(0.5, 0), coordsA=ax_i.transAxes,
        xyB=(image_time_numeric, bx.get_ylim()[1]), coordsB=bx.transData, 
        ls='--')
    ax_i.add_artist(line)
    bx.axvline(image_time_numeric, c='k', ls='--', alpha=1)

    ax_i.annotate("Pulsating\naurora", 
            xy=(-136, 61.24), xytext=(-133.97, 60.27),
            arrowprops=dict(arrowstyle="->", color='yellow'), color='yellow')

microburst_arrow_times = matplotlib.dates.date2num([
    datetime(2007, 2, 14, 13, 30, 39),
    datetime(2007, 2, 14, 13, 30, 36),
    datetime(2007, 2, 14, 13, 30, 42),
    ])
bx.annotate("Microbursts", 
    xy=(microburst_arrow_times[1], 705), 
    xytext=(microburst_arrow_times[0], 1500),
    arrowprops=dict(arrowstyle="->"), color='k', ha='center')
bx.annotate("", 
    xy=(microburst_arrow_times[2], 953), 
    xytext=(microburst_arrow_times[0], 1500),
    arrowprops=dict(arrowstyle="->"), color='k', ha='center')  

plt.suptitle(f'Example THEMIS ASI-SAMPEX Conjunction | Footprint {map_alt=} km', fontsize=15)
plt.subplots_adjust(hspace=0.10, wspace=0.01, top=0.93, bottom=0.2, left=0.09, right=0.95)
plt.show()
pass