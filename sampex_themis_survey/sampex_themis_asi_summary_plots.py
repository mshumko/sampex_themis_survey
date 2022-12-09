from datetime import datetime, date, timedelta
import dateutil.parser
import string
import pathlib

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
import aacgmv2

from sampex_themis_survey import config
from sampex_themis_survey.footprint import SAMPEX_footprint

class Summary:
    def __init__(self, conjunction_name, plot_pad_s=15, n_images=4, 
                map_alt=110) -> None:
        """
        Make summary plots of THEMIS ASI-SAMPEX conjunctions using a
        pre-computed conjunction list.

        Parameters
        ----------
        conjunction_name: str
            The filename of the conjunction csv file. The file must be in
            the config['data_dir'] directory.
        plot_pad_s: float
            Controls how much to expand the plot window, in seconds, beyond
            the "start" and "end" columns in the conjunction list.
        n_images: int
            The number of ASI images to plot in the first row.
        map_alt: int
            The footprint map altitude in kilometers.  
        """
        self.conjunction_name = conjunction_name
        self.plot_pad_s = plot_pad_s
        self.n_images = n_images
        self.map_alt = map_alt
        self.sampex_x_labels = {
            'L':'L_Shell', 'MLT':'MLT', 'Geo Lat':'GEO_Lat', 'Geo Lon':'GEO_Long'
            }
        self.plot_save_dir = config['plots_dir'] / datetime.now().strftime('%Y%m%d')
        self._load_conjunctions()
        pass

    def loop(self):
        """
        Loop over every conjunction and make a summary plot.
        """
        self.current_date = date.min
        if not self.plot_save_dir.exists():
            self.plot_save_dir.mkdir(parents=True)
            print(f'Created a {str(self.plot_save_dir)} plot directory.')

        self._init_plot()

        for (i, self.row) in self.conjunction_list.iterrows():
            if self.current_date != self.row['start'].date:
                # Load the SAMPEX-HILT data for this day.
                self.hilt = sampex.HILT(self.row['start']).load()
                self.current_date = self.row['start'].date
            print(f'Processing {self.row["start"]} ({i}/{self.conjunction_list.shape[0]})')

            self.time_range = [
                self.row['start'] - timedelta(seconds=self.plot_pad_s),
                self.row['end'] + timedelta(seconds=self.plot_pad_s)
            ]
            dt = (self.time_range[1]-self.time_range[0]).total_seconds()/(self.n_images-1)
            image_times = [self.time_range[0] + j*timedelta(seconds=dt) 
                for j in range(self.n_images)]

            self._get_footprint()
            # actual because self.image_times are calculates regardless of if there is an
            # image or not.
            self.actual_image_times = self._plot_images(image_times)
            self._plot_footprint()
            self._plot_hilt()
            self._format_subplots()
            self._plot_connecting_lines()
            plt.suptitle(
                f'{self.row["start"].date()} | '
                f'THEMIS {self.row["asi"].upper()}-SAMPEX Conjunction | '
                f'{self.map_alt} km footprint', 
                fontsize=15)

            plot_save_name = (
                f'{self.row["start"].strftime("%Y%m%d")}_'
                f'{self.row["start"].strftime("%H%M%S")}_'
                f'{self.row["end"].strftime("%H%M%S")}_'
                f'themis_{self.row["asi"].lower()}_'
                f'sampex_conjunction.png'
                )
            plt.savefig(self.plot_save_dir / plot_save_name)
            with open(self.last_movie_path, 'w') as f:
                f.write(self.row['start'].isoformat())
            self._clear_plot()  # After every iteration.
        return

    def _init_plot(self):
        """
        Initialize a plot with two rows. First row for self.n_images number of
        images and the second row the SAMPEX HILT data.
        """
        self.fig = plt.figure(figsize=(10, 5.5))
        spec = gridspec.GridSpec(
            nrows=3, ncols=self.n_images, figure=self.fig, height_ratios=(2, 1, 1)
            )
        self.ax = self.n_images*[None]
        for i, ax_i in enumerate(self.ax):
            # Using direct self.ax[i] reference so the output is 
            # saved into the self.ax list.
            self.ax[i] = self.fig.add_subplot(spec[0, i])
            self.ax[i].get_xaxis().set_visible(False)
            self.ax[i].get_yaxis().set_visible(False)

        self.bx = self.fig.add_subplot(spec[1, :])
        self.bx.get_xaxis().set_visible(False)
        self.cx = self.fig.add_subplot(spec[-1, :], sharex=self.bx)
        plt.subplots_adjust(
            hspace=0.10, wspace=0.01, top=0.93, bottom=0.2, left=0.09, right=0.95
            )
        return

    def _format_subplots(self):
        """
        Format the bx subplot with x-axis labels and connect to the x-axis
        formatter function.
        """
        self.bx.set_ylabel(f'mlat [deg]')
        self.bx.text(
            0, 0.99, f'({string.ascii_uppercase[self.n_images]}) '
            f'THEMIS-{self.row["asi"].upper()} keogram along footprint', 
            va='top', transform=self.bx.transAxes, weight='bold', fontsize=12,
            color='white'
            )
            
        self.cx.xaxis.set_major_formatter(FuncFormatter(self._format_fn))
        self.cx.xaxis.set_major_locator(matplotlib.dates.SecondLocator(interval=10))
        self.cx.set_xlabel('\n'.join(['Time'] + list(self.sampex_x_labels.keys())))
        self.cx.xaxis.set_label_coords(-0.07, -0.09)
        self.cx.text(
            0, 0.99, f'({string.ascii_uppercase[self.n_images+1]}) SAMPEX-HILT', 
            va='top', transform=self.cx.transAxes, weight='bold', fontsize=12
            )
        self.cx.set_xlim(self.actual_image_times[0], self.actual_image_times[-1])
        self.cx.set_ylabel(f'>1 MeV electrons\n[counts/20 ms]')
        return

    def _clear_plot(self):
        """
        Clears the plot for the next iteration.
        """
        for ax_i in self.ax:
            ax_i.clear()
        self.bx.clear()
        self.cx.clear()
        return

    def _plot_images(self, image_times):
        nearest_asi_image_times = []

        # Calculate the lat and lon bounds for pixels above 10 degree elevation.
        skymap = asilib.load_skymap('THEMIS', self.row['asi'], self.row['start'])
        # assert (
        #     self.map_alt in skymap['FULL_MAP_ALTITUDE'] / 1000
        # ), f'{self.map_alt} km is not in skymap calibration altitudes: {skymap["FULL_MAP_ALTITUDE"]/1000} km'
        # alt_index = np.where(skymap['FULL_MAP_ALTITUDE'] / 1000 == self.map_alt)[0][0]
        # idx_horizon = np.where(
        #     (skymap['FULL_ELEVATION'] < 50) & np.isnan(skymap['FULL_ELEVATION'])
        #     )
        # lat_map = skymap['FULL_MAP_LATITUDE'][alt_index, :, :].copy()
        # lon_map = skymap['FULL_MAP_LONGITUDE'][alt_index, :, :].copy()
        # lat_map[idx_horizon] = np.nan
        # lon_map[idx_horizon] = np.nan
        # lat_bounds = [np.nanmin(lat_map), np.nanmax(lat_map)]
        # lon_bounds = [np.nanmin(lon_map), np.nanmax(lon_map)]

        asilib.plot_keogram('THEMIS', self.row['asi'], self.time_range, 
            ax=self.bx, map_alt=self.map_alt, title=False, aacgm=True,
            path=self.footprint.loc[:, ['GEO_Lat', 'GEO_Long']].to_numpy()
            )

        lat_bounds = (
            skymap['SITE_MAP_LATITUDE']-6,
            skymap['SITE_MAP_LATITUDE']+6
            )
        lon_bounds = (
            skymap['SITE_MAP_LONGITUDE']-12,
            skymap['SITE_MAP_LONGITUDE']+12
            )

        z = zip(self.ax, image_times, string.ascii_uppercase[:self.n_images])
        for i, (ax_i, image_time, subplot_letter) in enumerate(z):
            asilib.make_map(lat_bounds=lat_bounds, lon_bounds=lon_bounds, ax=ax_i)
            try:
                t, _, _, _, _ = asilib.plot_map('THEMIS', self.row['asi'], image_time, self.map_alt, 
                    ax=ax_i, asi_label=False, color_bounds=None, pcolormesh_kwargs={'rasterized':True})
            except AssertionError as err:
                if '0 number of time stamps were found within 3 seconds' in str(err):
                    continue
                else:
                    raise
            nearest_asi_image_times.append(t)
            ax_i.text(
                0, 0.98, f'({subplot_letter}) {t.strftime("%H:%M:%S")}',
                transform=ax_i.transAxes, weight='bold', fontsize=12, va='top',
                color='purple'
                )
        return nearest_asi_image_times
    
    def _plot_footprint(self):
        # On the images
        for (ax_i, t) in zip(self.ax, self.actual_image_times):
            ax_i.plot(self.footprint.loc[:, 'GEO_Long'], self.footprint.loc[:, 'GEO_Lat'], 'r:')

            dt = [(t - t_fp).total_seconds() for t_fp in self.footprint.index]
            footprint_time = self.footprint.index[np.argmin(np.abs(dt))]
            ax_i.scatter(
                self.footprint.loc[footprint_time, 'GEO_Long'], 
                self.footprint.loc[footprint_time, 'GEO_Lat'],
                c='red', s=150, marker='.',
                )

        # On the keogram
        self.footprint['mlat'] = aacgmv2.convert_latlon_arr(
            self.footprint.loc[:, 'GEO_Lat'],
            self.footprint.loc[:, 'GEO_Long'],
            self.map_alt,
            self.time_range[0],
            method_code="G2A",
        )[0]
        self.bx.plot(self.footprint.index, self.footprint.mlat, 'r:')
        return

    def _plot_hilt(self):
        hilt_flt = self.hilt.loc[
            (self.hilt.index > self.time_range[0]) &
            (self.hilt.index <= self.time_range[1])
            ]
        self.cx.plot(hilt_flt.index, hilt_flt['counts'], c='k')
        self.cx.xaxis.set_minor_locator(matplotlib.dates.SecondLocator())
        return

    def _get_footprint(self):
        # Now load, map, and interpolate the SAMPEX attitude data
        asi_times, _ = asilib.load_image('THEMIS', self.row['asi'], 
            # time_range[0] is padded to include the first SAMPEX attitude point.
            time_range=(
                self.time_range[0]-pd.Timedelta(seconds=3), self.time_range[1]
                )
            )
        asi_times_df = pd.DataFrame(index=asi_times)
        self.footprint = SAMPEX_footprint(self.time_range[0]).map(alt=self.map_alt)

        # The merge will double the time cadence to 3 seconds and the new times will have
        # an associated NaNs. df.interpolate will interpolate over the NaN values.
        self.footprint = pd.merge_asof(
            asi_times_df, self.footprint, 
            left_index=True, right_index=True, 
            tolerance=pd.Timedelta(seconds=2), 
            direction='nearest'
            )
        self.footprint = self.footprint.interpolate(
            method='time', limit_area='inside'
            ).dropna()
        return

    def _plot_connecting_lines(self):
        """ 
        Draw lines connecting the subplots in the 0th row with the subplot 
        in the 1st row. 
        """
        z = zip(self.ax, matplotlib.dates.date2num(self.actual_image_times))
        for ax_i, image_time_numeric in z:
            line = matplotlib.patches.ConnectionPatch(
                xyA=(0.5, 0), coordsA=ax_i.transAxes,
                xyB=(image_time_numeric, self.bx.get_ylim()[1]), coordsB=self.bx.transData, 
                ls='--')
            ax_i.add_artist(line)
            self.bx.axvline(image_time_numeric, c='k', ls='--', alpha=1) 
            self.cx.axvline(image_time_numeric, c='k', ls='--', alpha=1) 
        return

    def _format_fn(self, tick_val, tick_pos):
        """
        The tick magic happens here. pyplot gives it a tick time, and this function 
        returns the closest label to that time. Read docs for FuncFormatter().
        """
        # Find the nearest time within 6 seconds (the cadence of the SAMPEX attitude files)
        tick_time = matplotlib.dates.num2date(tick_val).replace(tzinfo=None)
        i_min_time = np.argmin(np.abs(self.footprint.index - tick_time))
        if np.abs(self.footprint.index[i_min_time] - tick_time).total_seconds() > 3:
            return tick_time.strftime('%H:%M:%S')
            # raise ValueError(f'Nearest timestamp to tick_time is more than 6 seconds away')
        pd_index = self.footprint.index[i_min_time]
        # Cast np.array as strings so that it can insert the time string. 
        values = self.footprint.loc[pd_index, self.sampex_x_labels.values()].to_numpy().round(2).astype(str)
        values = np.insert(values, 0, tick_time.strftime('%H:%M:%S'))
        label = '\n'.join(values)
        return label

    def _load_conjunctions(self):
        """
        Load conjunction csv dataset. If last_summary_movie.txt exists, filter 
        the conjunctions by date and time to only events after the date in 
        last_summary_movie.txt.
        """
        conjunction_path = config['data_dir'] / self.conjunction_name

        self.conjunction_list = pd.read_csv(conjunction_path)
        self.conjunction_list['start'] = pd.to_datetime(self.conjunction_list['start'])
        self.conjunction_list['end'] = pd.to_datetime(self.conjunction_list['end'])

        self.last_movie_path = self.plot_save_dir / 'last_summary_movie.txt'
        if self.last_movie_path.exists():
            with open(self.last_movie_path, 'r') as f:
                last_movie_time = dateutil.parser.parse(f.read())
            self.conjunction_list = self.conjunction_list[
                self.conjunction_list['start'] > last_movie_time
                ]
        return

if __name__ == '__main__':
    conjunction_name = f'sampex_themis_asi_conjunctions_filtered.csv'
    s = Summary(conjunction_name)
    s.loop()