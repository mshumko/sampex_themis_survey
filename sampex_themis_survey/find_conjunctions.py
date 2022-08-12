"""
Find the conjunctions between SAMPEX and THEMIS ASIs
"""
from datetime import datetime
import pathlib

import pandas as pd
import numpy as np
import asilib

from sampex_survey import config
from sampex_survey.footprint import SAMPEX_footprint 
from proton_aurora.load import sampex

themis_imagers = asilib.themis_info()
themis_url = 'https://data.phys.ucalgary.ca/sort_by_project/THEMIS/asi/stream0/'


# Prepare the data/conjunction directory.
save_dir = pathlib.Path(config['data_dir'], 'data', 'conjunctions')
if not save_dir.exists():
    save_dir.mkdir(parents=True)
    print(f'Made {save_dir} directory.')
else:
    # Remove the csv files in the conjunction folder to avoid duplicate rows.
    for file in save_dir.glob('*.csv'):
        if file.is_file():
            file.unlink()
    print('Cleaned up the conjunction folder.')


days = pd.date_range(start='2006-12-01', end='2013-01-01', freq='D')
for day in days:
    # Look for a HILT file and skip the day if none exists.
    try:
        hilt = sampex.Load_HILT(day)
        # hilt.resolve_counts_state4()
    except (AssertionError, ValueError) as err:
        if '0 matched HILT files found.' in str(err):
            continue
        elif 'The SAMPEX HILT data is not in order' in str(err):
            continue
        else:
            raise
    
    print(f'Processing SAMPEX-THEMIS ASI conjunctions on {day.date()}')
    # Calculate the SAMPEX footprints in the LLA coordinates.
    footprint = SAMPEX_footprint(day)
    footprint.map_footprint()

    # What ASI was below SAMPEX?
    for location_code in themis_imagers['location_code']:
        try:
            imager = asilib.themis(location_code, time=day, load_images=False)
        except (Exception, AssertionError) as err:
            if 'Invalid SIGNATURE' in str(err):  # Poorly formatted save file.
                continue
            elif 'Only one href is allowed but' in str(err):
                continue
            else:
                raise
        c2 = asilib.Conjunction(
            imager, 
            footprint.attitude.index, 
            footprint.attitude.loc[:, ['GEO_Lat', 'GEO_Long', 'Altitude']].to_numpy())
        conjunction_df = c2.find()
        conjunction_df['hilt_data'] = False
        conjunction_df['asi_data'] = False

        for index, row in conjunction_df.iterrows():
            # Was there HILT data during this conjunction?
            idx = np.where(
                (hilt.hilt.index > row['start']) &
                (hilt.hilt.index < row['end'])
            )[0]
            if len(idx):
                conjunction_df.loc[index, 'hilt_data'] = True

            # Was there ASI data during this conjunction?
            download = asilib.io.download.Downloader(themis_url)
            url_subdirectories = [
                str(day.year), 
                str(day.month).zfill(2), 
                str(day.day).zfill(2), 
                f'{location_code.lower()}*', 
                f'ut{str(row["start"].hour).zfill(2)}', 
                ]
            filename = row["start"].strftime(f'%Y%m%d_%H%M_{location_code.lower()}*.pgm.gz')
            try:
                asi_url = download.find_url(subdirectories=url_subdirectories, filename=filename)
            except (FileNotFoundError, AssertionError) as err:
                if 'does not contain any hyper references containing' in str(err):
                    continue
                elif 'Only one href is allowed but' in str(err):
                    continue
                else:
                    raise
            if len(asi_url):
                conjunction_df.loc[index, 'asi_data'] = True
            
            # Save to file.
            save_name = f'sampex_themis_{location_code.lower()}_conjunctions.csv'
            save_path = save_dir / save_name

            if save_path.exists():
                conjunction_df.to_csv(save_path, mode='a', header=False, index=False)
            else:
                conjunction_df.to_csv(save_path, index=False)


# Finally merge the conjunction files into one.
file_paths = save_dir.rglob(f'sampex_themis_*_conjunctions.csv')
merged_conjunctions = pd.DataFrame(columns=['start', 'end', 'hilt_data', 'asi_data', 'asi'])

for file_path in file_paths:
    df = pd.read_csv(file_path)
    df['asi'] = file_path.name.split('_')[2]

    merged_conjunctions = pd.concat([merged_conjunctions, df])

# Sort the conjunctions by time.
merged_conjunctions['start'] = pd.to_datetime(merged_conjunctions['start'])
merged_conjunctions = merged_conjunctions.sort_values('start')
merged_conjunctions.to_csv(save_dir.parents[0] / f'sampex_themis_asi_conjunctions.csv', 
    index=False)

filtered_conjunctions = merged_conjunctions[
    (merged_conjunctions['asi_data'] == True) & (merged_conjunctions['hilt_data'] == True)
    ]
filtered_conjunctions.to_csv(save_dir.parents[0] / f'sampex_themis_asi_conjunctions_filtered.csv', 
    index=False)