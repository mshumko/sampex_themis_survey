"""
Make movies of ELFIN-THEMIS ASI conjunctions.
"""
import pathlib
from datetime import datetime, date, timedelta
import dateutil.parser

import pandas as pd
import numpy as np
import asilib
import matplotlib.pyplot as plt
import matplotlib.dates
from matplotlib.ticker import FuncFormatter
import sampex

from sampex_survey import config
from sampex_survey.footprint import SAMPEX_footprint 


alt = 110  # km
box = (10, 10)  # km
x_labels = {'L':'L_Shell', 'MLT':'MLT', 'Geo Lat':'GEO_Lat', 'Geo Lon':'GEO_Long'}

conjunction_dir = pathlib.Path(config.PROJECT_DIR, 'data')
conjunction_path = conjunction_dir / f'sampex_themis_asi_themis_aurorax_conjunctions.xlsx'

conjunction_list = pd.read_excel(conjunction_path, skiprows=1)
conjunction_list['start'] = pd.to_datetime(conjunction_list['Start Time (UTC)'])
conjunction_list['end'] = pd.to_datetime(conjunction_list['End Time (UTC)'])

current_date = date.min

last_movie_path = pathlib.Path(conjunction_dir, 'last_aurorax_summary_movie.txt')
if last_movie_path.exists():
    with open(last_movie_path, 'r') as f:
        last_movie_time = dateutil.parser.parse(f.read())
    conjunction_list = conjunction_list[conjunction_list['start'] > last_movie_time]

def format_fn(tick_val, tick_pos):
    """
    The tick magic happens here. pyplot gives it a tick time, and this function 
    returns the closest label to that time. Read docs for FuncFormatter().
    """
    # Find the nearest time within 6 seconds (the cadence of the SAMPEX attitude files)
    tick_time = matplotlib.dates.num2date(tick_val).replace(tzinfo=None)
    i_min_time = np.argmin(np.abs(footprint.attitude.index - tick_time))
    if np.abs(footprint.attitude.index[i_min_time] - tick_time).total_seconds() > 6:
        return tick_time.strftime('%H:%M:%S')
    pd_index = footprint.attitude.index[i_min_time]
    # Cast np.array as strings so that it can insert the time string. 
    values = footprint.attitude.loc[pd_index, x_labels.values()].to_numpy().round(2).astype(str)
    values = np.insert(values, 0, pd_index.strftime('%H:%M:%S'))
    label = '\n'.join(values)
    return label

for _, row in conjunction_list.iterrows():
    print(f'Processing {row["start"]}')
    if current_date != row['start'].date:
        try:
            hilt = sampex.HILT(row['start']).load()
        except FileNotFoundError as err:
            if 'does not contain any hyper references' in str(err):
                continue
            raise
        current_date = row['start'].date

    asi_location = row['Conjunction Between'].split('and')[0].rstrip().split(' ')[1]
    time_range = [
        row['start'] - timedelta(minutes=0.25),
        row['end'] + timedelta(minutes=0.25)
    ]
    # Create an Imager object
    img = asilib.themis(asi_location, time_range=time_range, alt=alt)
    # Load, filter, and map the SAMPEX footprint
    footprint = SAMPEX_footprint(row['start'])
    footprint.attitude = footprint.attitude.loc[
        (footprint.attitude.index >= time_range[0]) &
        (footprint.attitude.index <= time_range[1]), :
    ]
    footprint.map(alt=alt)

    c = asilib.Conjunction(
        img, 
        footprint.attitude.index, 
        footprint.attitude.loc[:, ['GEO_Lat', 'GEO_Long', 'Altitude']].to_numpy()
        )
    try:
        c.resample()
    except ValueError as err:
        if 'No imager time stamps to interpolate over.' in str(err):
            continue
        else:
            raise
    sat_azel, asi_pixels = c.map_lla_azel()
    equal_area_gen = c.equal_area_gen(box=box)
    mask_gen = c.equal_area_gen(box=(10, 10))

    mean_asi_intensity = -np.ones(c.sat.shape[0])
    for i, ((_, image), mask) in enumerate(zip(img, equal_area_gen)):
        mean_asi_intensity[i] = np.nanmean(image*mask)
    
    fig, ax = plt.subplots(3, gridspec_kw={'height_ratios':[3, 1, 1]}, figsize=(6, 10))
    ax[1].sharex(ax[2])  # Connect the two subplots to remove the extra time axis.

    hilt_copy = hilt[(hilt.index >= time_range[0]) & (hilt.index <= time_range[1])]

    img_gen = img.animate_fisheye_gen(ax=ax[0], overwrite=True)
    mask_gen = c.equal_area_gen(box=box)

    for i, ((image_time, image, _, im), mask) in enumerate(zip(img_gen, mask_gen)):
        ax[1].clear()
        ax[2].clear()
        ax[1].xaxis.set_visible([])

        ax[0].plot(asi_pixels[:, 0], asi_pixels[:, 1], 'r:')
        ax[0].scatter(asi_pixels[i, 0], asi_pixels[i, 1], s=10, c='r')

        # Plot the equal area
        mask[np.isnan(mask)] = 0  # To play nice with plt.contour()
        ax[0].contour(mask, levels=[0.99], colors=['yellow'])

        ax[1].plot(c.sat.index, mean_asi_intensity, color='k')
        ax[2].plot(hilt_copy.index, hilt_copy.counts, color='k')
        ax[1].axvline(image_time, c='k')
        ax[2].axvline(image_time, c='k')
        ax[1].set_ylabel(f'Mean ASI intensity\n{box} km area')
        ax[2].set_ylabel(f'HILT')
        ax[2].set_xlabel(f'Time')
        ax[2].set_xlim(*time_range)

        ax[-1].xaxis.set_major_formatter(FuncFormatter(format_fn))
        ax[-1].set_xlabel('\n'.join(['Time'] + list(x_labels.keys())))
        ax[-1].xaxis.set_minor_locator(matplotlib.dates.SecondLocator())
        ax[-1].xaxis.set_label_coords(-0.1, -0.06)

        # ax[2].set_yscale('log')
        plt.suptitle(
            f"Conjunction between SAMPEX & THEMIS-{asi_location.upper()}\n"
            f'{img._data["time_range"][0].strftime("%Y-%m-%d")} | '
            f'{alt} km footprint'
            )
        # plt.tight_layout()
        plt.subplots_adjust(hspace=0.01, wspace=0.01, top=0.95, bottom=0.10, left=0.15, right=0.95)


    plt.close()
    with open(last_movie_path, 'w') as f:
        f.write(row['start'].isoformat())