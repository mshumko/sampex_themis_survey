"""
Plot the instantaneous locations of the THEMIS ASIs, and THEMIS probe footprints 
during  the 2008-03-04 triple conjunction.
"""
from datetime import datetime, date, timedelta
from importlib.resources import path
import string
import dateutil.parser
import pathlib
import logging

import numpy as np
import pandas as pd
import asilib
import sampex
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from matplotlib.ticker import FuncFormatter
import matplotlib.colors
import matplotlib.dates

import sampex_themis_survey
from sampex_themis_survey.themis.footprint import THEMIS_footprint
from sampex_themis_survey.footprint import SAMPEX_footprint
from sampex_themis_survey.themis.sst import SST
from sampex_themis_survey.themis.fbk import FBK


class Themis_Themis_ASI:
    """
    Make summary plots of conjunctions between the THEMIS probes and THEMIS ASIs.
    """
    def __init__(self, c_filename, time_window_sec=3600/2, n_images=3, map_alt=110) -> None:
        self.c_filename = c_filename
        self.time_window_sec = time_window_sec
        self.n_images = n_images
        self.map_alt = map_alt

        self.save_dir = pathlib.Path(sampex_themis_survey.config['code_dir'], 'plots', 
            c_filename.split('.')[0])
        self.save_dir.mkdir(parents=True, exist_ok=True)

        self._load()
        return

    def _load(self):
        self.c_path = pathlib.Path(sampex_themis_survey.config['code_dir'], 'data', self.c_filename)
        self.c_df = pd.read_excel(self.c_path, skiprows=1)
        self.c_df['start'] = pd.to_datetime(self.c_df['Start Time (UTC)'])
        self.c_df['end'] = pd.to_datetime(self.c_df['End Time (UTC)'])

        # Don't reanalyze conjunctions if last_summary_plot.txt exists.
        self.last_plot_date_path = pathlib.Path(self.save_dir, 'last_summary_plot.txt')
        if self.last_plot_date_path.exists():
            with open(self.last_plot_date_path, 'r') as f:
                last_plot_time = dateutil.parser.parse(f.read())
            self.c_df = self.c_df[self.c_df['start'] > last_plot_time]
        return

    def loop(self):
        """
        Loop over every conjunction and load the data.
        """
        self._create_empty_subplots()
        current_date = date.min
        for _, self.row in self.c_df.iterrows():
            print(f'Processing {self.row["start"]}')
            if current_date != self.row['start'].date:
                try:
                    self.hilt = sampex.HILT(self.row['start']).load()
                except FileNotFoundError as err:
                    if 'does not contain any hyper references' in str(err):
                        continue
                    raise
                current_date = self.row['start'].date

            self.asi_location = self.row['Conjunction Between'].split('and')[0].rstrip().split(' ')[1]
            self.sc_id = self.row['Conjunction Between'].split('and')[1].rstrip().split('-')[1]
            self.time_range = [
                self.row['start'] - timedelta(seconds=self.time_window_sec/2),
                self.row['end'] + timedelta(seconds=self.time_window_sec/2)
            ]
            

            self.plot()

            with open(self.last_plot_date_path, 'w') as f:
                f.write(self.row['start'].isoformat())
        
        return

    def plot(self, save=True):
        """
        Creates the subplot layout and dispatches the plotting to the other methods.
        """
        for ax_i in self.ax:
            ax_i.clear()
        for ax_i in [self.bx, self.cx, self.dx, self.ex]:
            ax_i.clear()

        # -1 so that image_time[-1] == time_range[1]
        dt = (self.time_range[1]-self.time_range[0])/(self.n_images-1)
        self.image_times = [self.time_range[0] +  i*dt for i in range(self.n_images)]
        self._plot_asi_images(self.ax)
        try:
            self._plot_themis_footprint(self.ax)
        except AssertionError as err:
            if 'The time stamp difference' in str(err):
                return
            else:
                raise
        self._plot_keogram(self.bx)

        # The try-except blocks in case one of the instruments did not have data
        # on that day.
        try:
            self._plot_sst_e(self.cx)
        except FileNotFoundError as err:
            if 'does not contain any hyper references containing' in str(err):
                return
            else:
                raise
        try:
            self._plot_sst_p(self.dx)
        except FileNotFoundError as err:
            if 'does not contain any hyper references containing' in str(err):
                return
            else:
                raise
        try:
            self._plot_fbk(self.ex)
        except FileNotFoundError as err:
            if 'does not contain any hyper references containing' in str(err):
                return
            else:
                raise

        self.ex.set_xlabel('Time [HH:MM]')
        time_format = matplotlib.dates.DateFormatter('%H:%M')
        self.ex.xaxis.set_major_formatter(time_format)

        # Annotate and clean up the plot.
        plot_labels = (
            'THEMIS-ASI keogram',
            'THEMIS-SST electrons',
            'THEMIS-SST ions',
            'THEMIS-FBK'
        )
        subplots = [self.bx, self.cx, self.dx, self.ex]
        label_colors = ['w', 'k', 'k', 'k']
        z = zip(subplots, plot_labels, label_colors)
        for i, (ax_i, plot_label, label_color) in enumerate(z, start=self.n_images):
            ax_i.text(0, 0, f'({string.ascii_uppercase[i]}) {plot_label}', 
                transform=ax_i.transAxes, va='bottom', color=label_color, fontsize=15)

        plt.suptitle(f'Example Triple Conjunction | THEMIS Probe-THEMIS ASI', fontsize=15)

        # plt.legend()
        plt.subplots_adjust(wspace=0.02, hspace=0.07, left=0.055, right=0.92, top=0.943)
        if save:
            save_time = self.row['start'].strftime("%Y%m%d_%H%M%S")
            filename = f'{save_time}_themis_probe_{self.sc_id}_themis_asi_{self.asi_location}_conjunction.png'
            plt.savefig(self.save_dir / filename, dpi=300)
        else:
            plt.show()
        return

    def _create_empty_subplots(self):
        """
        Create the empty subplot layout and make the subplots
        as class attributes. 
        """
        self.fig = plt.figure(figsize=(12, 9))
        self.spec = gridspec.GridSpec(
            nrows=5, 
            ncols=self.n_images, figure=self.fig, height_ratios=(2, 1, 1, 1, 1)
            )
        # ASI map subplots
        self.ax = self.n_images*[None]
        for i in range(self.n_images):
            self.ax[i] = self.fig.add_subplot(self.spec[0, i])
            self.ax[i].get_xaxis().set_visible(False)
            self.ax[i].get_yaxis().set_visible(False)
        # Keogram
        self.bx = self.fig.add_subplot(self.spec[1, :])
        self.bx.tick_params(axis="x", labelbottom=False)
        # THEMIS-SST electrons
        self.cx = self.fig.add_subplot(self.spec[2, :], sharex=self.bx) 
        self.cx.tick_params(axis="x", labelbottom=False)
        # THEMIS-SST protons
        self.dx = self.fig.add_subplot(self.spec[3, :], sharex=self.bx) 
        self.dx.tick_params(axis="x", labelbottom=False) 
        # THEMIS-FBK magnetic spectrum. No need for the tick_params
        # command since it is the bottom subplot.
        self.ex = self.fig.add_subplot(self.spec[4, :], sharex=self.bx)
        return


    def _plot_asi_images(self, _ax):
        """
        Plot a sequence of ASI images projected onto a geographic map in the first row
        of the subplot grid.
        """
        z = zip(_ax, self.image_times, string.ascii_uppercase[:self.n_images])

        skymap = asilib.load_skymap('THEMIS', self.asi_location, self.image_times[0])
        self.lat_bounds = (
            skymap['SITE_MAP_LATITUDE']-6,
            skymap['SITE_MAP_LATITUDE']+6
            )
        self.lon_bounds = (
            skymap['SITE_MAP_LONGITUDE']-12,
            skymap['SITE_MAP_LONGITUDE']+12
            )

        for i, (ax_i, image_time, subplot_letter) in enumerate(z):
            asilib.make_map(lat_bounds=self.lat_bounds, lon_bounds=self.lon_bounds,
                ax=ax_i, land_color='grey')

            try:
                # Update image_times with the actual image timestamp.
                self.image_times[i], _, _, _, _ = asilib.plot_map(
                    'THEMIS', self.asi_location, image_time, self.map_alt, 
                    ax=ax_i, asi_label=False, color_bounds=None, pcolormesh_kwargs={'rasterized':True},
                    time_thresh_s=3
                    )
            except AssertionError as err:
                if '0 number of time stamps were found within' in str(err):
                    continue
                else:
                    raise
            
            ax_i.text(0, 0, f'({subplot_letter}) {image_time.strftime("%H:%M:%S")}', 
                transform=ax_i.transAxes, va='bottom', color='k', fontsize=15)
        return

    def _plot_themis_footprint(self, _ax):
        """
        Plot the THEMIS probe footprint location in the first row of the subplot grid.
        """
        f = THEMIS_footprint(self.sc_id, self.time_range)
        f.map(alt=self.map_alt, model='T89', mag_input={'Kp':1.667})
        
        for time, ax_i in zip(self.image_times, _ax):
            dt = np.abs([(time-ti).total_seconds() for ti in f.time])
            idt = np.argmin(dt)
            assert (time-f.time[idt]).total_seconds() <= 60, (
                f'The time stamp difference, {(time-f.time[idt]).total_seconds()} '
                f's, is too large.')
            ax_i.scatter(f.footprint[idt, 1], f.footprint[idt, 0], 
                label=f'THEMIS-{self.sc_id}',
                color='r', s=100, marker='x')
        return

    def _plot_keogram(self,_ax):
        """
        Plots a keogram on the 2nd row of the subplot grid.
        """
        try:
            asilib.plot_keogram('THEMIS', self.asi_location, self.time_range, 
                ax=_ax, map_alt=self.map_alt, aacgm=True)
        except ValueError as err:
            if 'time stamps are out of order' in str(err):
                print(err)
                return
            else:
                raise
        _ax.set_xlim(self.time_range)
        _ax.set_title('')
        _ax.set_ylabel('Magnetic lat [deg]')
        # self.bx.set_ylim(self.lat_bounds)
        return

    def _plot_sst_e(self, _ax):
        """ 
        Plot the THEMIS-SST electrons.
        """ 
        # Electrons
        ax_c = _ax.inset_axes([1.01, 0, 0.02, 1], transform=_ax.transAxes)
        s = SST(self.sc_id, self.time_range[0], species='e')
        s.load()
        _, p = s.spectrum(ax=_ax, 
            pcolormesh_kwargs={'norm':matplotlib.colors.LogNorm(vmin=1E2, vmax=5E5)}
            )
        plt.colorbar(p, ax=_ax, cax=ax_c, label='Electron flux')
        _ax.set_xlim(self.time_range)
        _ax.set_ylabel('Energy [keV]')
        _ax.set_yscale('log')
        return

    def _plot_sst_p(self, _ax):
        """ 
        Plot the THEMIS-SST protons.
        """
        ax_c = _ax.inset_axes([1.01, 0, 0.02, 1], transform=_ax.transAxes)
        s = SST(self.sc_id, self.time_range[0], species='i')
        s.load()
        _, p = s.spectrum(ax=_ax, 
            pcolormesh_kwargs={'norm':matplotlib.colors.LogNorm(vmin=0.1, vmax=1E6)}
            )
        plt.colorbar(p, cax=ax_c, label='ion flux')
        _ax.set_ylabel('Energy [keV]')
        _ax.set_yscale('log')
        return

    def _plot_fbk(self, _ax):
        """
        Plot the THEMIS-FBK wave data
        """
        ax_c = _ax.inset_axes([1.01, 0, 0.02, 1], transform=_ax.transAxes)
        s = FBK(self.sc_id, self.time_range[0])
        try:
            s.load()
        except ValueError as err:
            if 'No records found for variable' in str(err):
                return
            else:
                raise
        _, p = s.spectrum(ax=_ax, 
            pcolormesh_kwargs={'norm':matplotlib.colors.LogNorm(vmin=1E-4, vmax=0.1)})
        plt.colorbar(p, cax=ax_c, label=f'$nT^{{{2}}}/Hz$')
        _ax.set_ylabel('Frequency [Hz]') 
        _ax.set_yscale('log')
        _ax.set_ylim(1, 1E3)
        return


class Sampex_Themis_ASI(Themis_Themis_ASI):
    def __init__(self, c_filename, time_window_sec=120, n_images=3, map_alt=110) -> None:
        super().__init__(
            c_filename, time_window_sec=time_window_sec, n_images=n_images, map_alt=map_alt
            )
        self.x_labels = {'L':'L_Shell', 'MLT':'MLT', 'Geo Lat':'GEO_Lat', 'Geo Lon':'GEO_Long'}
        logging.basicConfig(filename='sampex_themis_asi_summary_plot.log', 
            encoding='utf-8', level=logging.INFO, format='%(asctime)s %(message)s')
        return

    def loop(self):
        """
        Loop over every conjunction and load the data.
        """
        self._create_empty_subplots()

        current_date = date.min
        for _, self.row in self.c_df.iterrows():
            self.asi_location = self.row['Conjunction Between'].split('and')[0].rstrip().split(' ')[1]
            self.sc_id = self.row['Conjunction Between'].split('and')[1].rstrip().split('-')[1]
            t0 = self.row['start'] + (self.row['end'] - self.row['start'])/2
            self.time_range = [
                t0 - timedelta(seconds=self.time_window_sec/2),
                t0 + timedelta(seconds=self.time_window_sec/2)
            ]
            print(f'Processing {self.row["start"]} conjunction between SAMPEX and THEMIS-{self.asi_location}.')

            if current_date != self.row['start'].date:
                try:
                    self._load_sampex()
                except FileNotFoundError as err:
                    if 'does not contain any hyper references' in str(err):
                        continue
                    raise
                current_date = self.row['start'].date

            self.plot()

            with open(self.last_plot_date_path, 'w') as f:
                f.write(self.row['start'].isoformat())
        
        return

    def plot_one_conjunction(self, asi_location, time_range):
        self.time_range = time_range
        self.asi_location = asi_location
        self._create_empty_subplots()

        self.row = {}
        self.row['start'] = self.time_range[0]
        self._load_sampex()
        self.plot(save=False)
        return

    def plot(self, save=True):
        """
        Creates the subplot layout and dispatches the plotting to the other methods.
        """
        for ax_i in self.ax:
            ax_i.clear()
        self.bx.clear()

        # -1 so that image_time[-1] == time_range[1]
        dt = (self.time_range[1]-self.time_range[0])/(self.n_images-1)
        self.image_times = [self.time_range[0] +  i*dt for i in range(self.n_images)]
        self._plot_asi_images(self.ax)
        self._plot_sampex_footprint(self.ax)

        try:
            self._plot_hilt(self.bx)
            sampex_loaded=True  # A hack to get around the 
        except (FileNotFoundError, KeyError) as err:
            sampex_loaded=False
            if 'does not contain any hyper references containing' in str(err):  # No file
                pass
            elif isinstance(err, KeyError):  # No HILT data throughout self.time_range.
                logging.info(f'SAMPEX-HILT did not have continuos data in {self.time_range=}')
                pass
            else:
                raise
        # try:
        #     self._plot_pet(self.cx)
        # except FileNotFoundError as err:
        #     if 'does not contain any hyper references containing' in str(err):
        #         return
        #     else:
        #         raise

        # Annotate and tidy up the plot.
        plot_labels = (
            'HILT',
        )
        subplots = [self.bx]
        z = zip(subplots, plot_labels)
        for i, (ax_i, plot_label) in enumerate(z, start=self.n_images):
            ax_i.text(0, 0.99, f'({string.ascii_uppercase[i]}) {plot_label}', 
                transform=ax_i.transAxes, va='top', color='k', fontsize=15)

        if sampex_loaded:
            self.bx.xaxis.set_major_formatter(FuncFormatter(self.format_xaxis))
            self.bx.xaxis.set_major_locator(matplotlib.dates.SecondLocator(interval=10))

        self.bx.set_xlabel('\n'.join(['Time'] + list(self.x_labels.keys())))
        self.bx.xaxis.set_label_coords(-0.07, -0.09)

        plt.suptitle(f'Example Triple Conjunction | THEMIS ASI-SAMPEX', fontsize=15)
        plt.subplots_adjust(hspace=0.10, wspace=0.01, top=0.93, bottom=0.2, left=0.09, right=0.95)
        if save:
            save_time = self.row['start'].strftime("%Y%m%d_%H%M%S")
            filename = f'{save_time}_sampex_themis_asi_{self.asi_location}_conjunction.png'
            plt.savefig(self.save_dir / filename, dpi=300)
        else:
            plt.show()
        return

    def _load_sampex(self):
        self.hilt = sampex.HILT(self.row['start']).load()

        # At a 6-second cadence
        self.footprint = SAMPEX_footprint(self.time_range[0]).map(alt=self.map_alt)
        # Resample to a higher cadence and linearly interpolate the footprint values for more 
        # accurate location. Necessary for 1) plotting the instaneous SAMPEX position at the imager time, 
        # and 2) the stacked x-axis labels.
        self.footprint = self.footprint.resample('1S').interpolate()
        return

    def _create_empty_subplots(self):
        """
        Create the empty subplot layout and make the subplots
        as class attributes. 
        """
        self.fig = plt.figure(figsize=(10, 5))
        self.spec = gridspec.GridSpec(nrows=2, ncols=self.n_images, figure=self.fig, height_ratios=(2, 1))

        # ASI map subplots
        self.ax = self.n_images*[None]
        for i in range(self.n_images):
            self.ax[i] = self.fig.add_subplot(self.spec[0, i])
            self.ax[i].get_xaxis().set_visible(False)
            self.ax[i].get_yaxis().set_visible(False)
        # HILT
        self.bx = self.fig.add_subplot(self.spec[1, :])
        return

    def _plot_asi_images(self, _ax):
        """
        Plot a sequence of ASI images projected onto a geographic map in the first row
        of the subplot grid.
        """
        z = zip(_ax, self.image_times, string.ascii_uppercase[:self.n_images])

        skymap = asilib.load_skymap('THEMIS', self.asi_location, self.image_times[0])
        self.lat_bounds = (
            skymap['SITE_MAP_LATITUDE']-6,
            skymap['SITE_MAP_LATITUDE']+6
            )
        self.lon_bounds = (
            skymap['SITE_MAP_LONGITUDE']-12,
            skymap['SITE_MAP_LONGITUDE']+12
            )

        for i, (ax_i, image_time, subplot_letter) in enumerate(z):
            asilib.make_map(lat_bounds=self.lat_bounds, lon_bounds=self.lon_bounds,
                ax=ax_i, land_color='grey')

            try:
                # Update image_times with the actual image timestamp.
                self.image_times[i], _, _, _, _ = asilib.plot_map(
                    'THEMIS', self.asi_location, image_time, self.map_alt, 
                    ax=ax_i, asi_label=False, color_bounds=None, pcolormesh_kwargs={'rasterized':True},
                    time_thresh_s=3
                    )
            except AssertionError as err:
                if '0 number of time stamps were found within' in str(err):
                    continue
                else:
                    raise
            
            ax_i.text(0, 0, f'({subplot_letter}) {image_time.strftime("%H:%M:%S")}', 
                transform=ax_i.transAxes, va='bottom', color='k', fontsize=15)
        return

    def _plot_sampex_footprint(self, _ax):
        for ax_i, t in zip(_ax, self.image_times):
            ax_i.plot(
                self.footprint.loc[self.time_range[0]:self.time_range[1], 'GEO_Long'], 
                self.footprint.loc[self.time_range[0]:self.time_range[1], 'GEO_Lat'], 
                'r:')

            nearest_time_i = np.argmin(np.abs(self.footprint.index - t))
            nearest_time = self.footprint.index[nearest_time_i]
            dt = np.abs(nearest_time - t).total_seconds()

            if dt > 3:
                # return ''
                raise ValueError(
                    f'Nearest timestamp to tick_time is {dt} seconds away '
                    '(more than the allowable 3-second threshold)'
                    )
            ax_i.scatter(
                self.footprint.loc[nearest_time, 'GEO_Long'], 
                self.footprint.loc[nearest_time, 'GEO_Lat'],
                c='red', s=150, marker='.',
                )
        return

    def _plot_hilt(self, _ax):
        filtered_hilt = self.hilt.loc[self.time_range[0]:self.time_range[1], :]
        _ax.plot(filtered_hilt.index, filtered_hilt['counts'], c='r')
        _ax.xaxis.set_minor_locator(matplotlib.dates.SecondLocator())
        _ax.set_xlim(*self.time_range)
        # _ax.set_yscale('log')
        _ax.set_ylabel(f'>1 MeV electrons\n[counts/20 ms]')
        return

    def _plot_pet(self, _ax):
        pet = sampex.PET(self.time_range[0]).load()
        # Factor of two since PET accumulated counts over 50 ms, but the timestamps are at
        # every 100 ms.
        pet['P1_Rate'] = 2*pet['P1_Rate']
        _ax.step(pet.index, pet['P1_Rate'], label='PET', where='post')
        _ax.set_yscale('log')
        _ax.set_ylabel(f'>400 keV electrons\n[counts/100 ms]')
        return

    def format_xaxis(self, tick_val, _):
        """
        The tick magic happens here. pyplot gives it a tick time, and this function 
        returns the closest label to that time. Read docs for FuncFormatter().
        """
        # Find the nearest time within 6 seconds (the cadence of the SAMPEX attitude files)
        tick_time = matplotlib.dates.num2date(tick_val).replace(tzinfo=None)
        i_min_time = np.argmin(np.abs(self.footprint.index - tick_time))
        dt = np.abs(self.footprint.index[i_min_time] - tick_time).total_seconds()
        if dt > 6:
            # return ''
            raise ValueError(
                    f'Nearest timestamp to tick_time is {dt} seconds away '
                    '(more than the allowable 6-second threshold)'
                    )
        pd_index = self.footprint.index[i_min_time]
        # Cast np.array as strings so that it can insert the time string. 
        values = self.footprint.loc[pd_index, self.x_labels.values()].to_numpy().round(2).astype(str)
        values = np.insert(values, 0, tick_time.strftime('%H:%M:%S'))
        label = '\n'.join(values)
        return label



if __name__ == '__main__':
    # filename = 'sampex_themis_asi_themis_aurorax_conjunctions_500_km.xlsx'
    # # s = Themis_Themis_ASI(filename)
    # # s.loop()
    # s = Sampex_Themis_ASI(filename)
    # # s.loop()
    # s.plot_one_conjunction('FSMI', 
    #     (datetime(2008, 3, 4, 5, 49, 0), datetime(2008, 3, 4, 5, 50, 40))
    #     )
    
    filename = 'sampex_themis_asi_rbsp_aurorax_conjunctions_1000_km.xlsx'
    s = Sampex_Themis_ASI(filename)
    s.loop()