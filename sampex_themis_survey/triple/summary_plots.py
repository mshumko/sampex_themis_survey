"""
Plot the instantaneous locations of the THEMIS ASIs, and THEMIS probe footprints 
during  the 2008-03-04 triple conjunction.
"""
from datetime import datetime, date, timedelta
from importlib.resources import path
import string
import dateutil.parser
import pathlib

import numpy as np
import pandas as pd
import asilib
import sampex
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import matplotlib.colors
import matplotlib.dates

import sampex_themis_survey
from sampex_themis_survey.themis.footprint import THEMIS_footprint
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
        # TODO: Create the subplots here and pass each subplot to the relevant method.
        self.fig = plt.figure(figsize=(12, 9))
        self.spec = gridspec.GridSpec(
            nrows=5, 
            ncols=self.n_images, figure=self.fig, height_ratios=(2, 1, 1, 1, 1)
            )
        # -1 so that image_time[-1] == time_range[1]
        dt = (self.time_range[1]-self.time_range[0])/(self.n_images-1)
        self.image_times = [self.time_range[0] +  i*dt for i in range(self.n_images)]
        self._plot_asi_images()
        try:
            self._plot_themis_footprint()
        except AssertionError as err:
            if 'The time stamp difference' in str(err):
                return
            else:
                raise
        self._plot_keogram()
        try:
            self._plot_sst_e()
            self._plot_sst_p()
            self._plot_fbk()
        except FileNotFoundError as err:
            if 'does not contain any hyper references containing' in str(err):
                return
            else:
                raise
        # self.spec.tight_layout(self.fig)
        # plt.show()
        plot_labels = (
            'THEMIS-SST electrons',
            'THEMIS-SST ions',
            'THEMIS-FBK'
        )
        subplots = [self.cx, self.dx, self.ex]
        z = zip(subplots, plot_labels)
        for i, (ax_i, plot_label) in enumerate(z, start=self.n_images+1):
            ax_i.text(0, 0, f'({string.ascii_uppercase[i]}) {plot_label}', 
                transform=ax_i.transAxes, va='bottom', color='k', fontsize=15)

        plt.suptitle(f'Example Triple Conjunction | THEMIS-THEMIS ASI-SAMPEX', fontsize=15)

        # plt.legend()
        plt.subplots_adjust(wspace=0.02, hspace=0.07, left=0.055, right=0.92, top=0.943)
        if save:
            save_time = self.row['start'].strftime("%Y%m%d_%H%M%S")
            filename = f'{save_time}_themis_probe_{self.sc_id}_themis_asi_{self.asi_location}_conjunction.png'
            plt.savefig(self.save_dir / filename, dpi=300)
        else:
            plt.show()
        return

    def _plot_asi_images(self):
        """
        Plot a sequence of ASI images projected onto a geographic map in the first row
        of the subplot grid.
        """
        self.ax = self.n_images*[None]
        z = zip(self.image_times, string.ascii_uppercase[:self.n_images])

        skymap = asilib.load_skymap('THEMIS', self.asi_location, self.image_times[0])
        self.lat_bounds = (
            skymap['SITE_MAP_LATITUDE']-6,
            skymap['SITE_MAP_LATITUDE']+6
            )
        self.lon_bounds = (
            skymap['SITE_MAP_LONGITUDE']-12,
            skymap['SITE_MAP_LONGITUDE']+12
            )

        for i, (image_time, subplot_letter) in enumerate(z):
            self.ax[i] = self.fig.add_subplot(self.spec[0, i])
            self.ax[i].get_xaxis().set_visible(False)
            self.ax[i].get_yaxis().set_visible(False)
            asilib.make_map(lat_bounds=self.lat_bounds, lon_bounds=self.lon_bounds,
                ax=self.ax[i], land_color='grey')

            try:
                # Update image_times with the actual image timestamp.
                self.image_times[i], _, _, _, _ = asilib.plot_map(
                    'THEMIS', self.asi_location, image_time, self.map_alt, 
                    ax=self.ax[i], asi_label=False, color_bounds=None, pcolormesh_kwargs={'rasterized':True},
                    time_thresh_s=3
                    )
            except AssertionError as err:
                if '0 number of time stamps were found within' in str(err):
                    continue
                else:
                    raise
            
            self.ax[i].text(0, 0, f'({subplot_letter}) {image_time.strftime("%H:%M:%S")}', 
                transform=self.ax[i].transAxes, va='bottom', color='k', fontsize=15)
        return

    def _plot_themis_footprint(self):
        """
        Plot the THEMIS probe footprint location in the first row of the subplot grid.
        """
        f = THEMIS_footprint(self.sc_id, self.time_range)
        f.map(alt=self.map_alt, model='T89', mag_input={'Kp':1.667})
        
        for time, ax_i in zip(self.image_times, self.ax):
            dt = np.abs([(time-ti).total_seconds() for ti in f.time])
            idt = np.argmin(dt)
            assert (time-f.time[idt]).total_seconds() <= 60, (
                f'The time stamp difference, {(time-f.time[idt]).total_seconds()} '
                f's, is too large.')
            ax_i.scatter(f.footprint[idt, 1], f.footprint[idt, 0], 
                label=f'THEMIS-{self.sc_id}',
                color='r', s=100, marker='x')
        return

    def _plot_keogram(self):
        """
        Plots a keogram on the 2nd row of the subplot grid.
        """
        self.bx = self.fig.add_subplot(self.spec[1, :])
        self.bx.tick_params(axis="x", labelbottom=False)
        try:
            asilib.plot_keogram('THEMIS', self.asi_location, self.time_range, 
                ax=self.bx, map_alt=self.map_alt, aacgm=True)
        except ValueError as err:
            if 'time stamps are out of order' in str(err):
                print(err)
                return
            else:
                raise
        self.bx.set_xlim(self.time_range)
        self.bx.set_title('')
        self.bx.text(0, 1, f'({string.ascii_uppercase[self.n_images+1]}) THEMIS ASI keogram', 
            transform=self.bx.transAxes, va='top', color='white', fontsize=15)
        self.bx.set_ylabel('Magnetic lat [deg]')
        # self.bx.set_ylim(self.lat_bounds)
        return

    def _plot_sst_e(self):
        """ 
        Plot the THEMIS-SST electrons.
        """
        self.cx = self.fig.add_subplot(self.spec[2, :], sharex=self.bx) 
        self.cx.tick_params(axis="x", labelbottom=False) 

        # Electrons
        cx_c = self.cx.inset_axes([1.01, 0, 0.02, 1], transform=self.cx.transAxes)
        s = SST(self.sc_id, self.time_range[0], species='e')
        s.load()
        _, p = s.spectrum(ax=self.cx, 
            pcolormesh_kwargs={'norm':matplotlib.colors.LogNorm(vmin=1E2, vmax=5E5)}
            )
        plt.colorbar(p, ax=self.cx, cax=cx_c, label='Electron flux')
        self.cx.set_xlim(self.time_range)
        self.cx.set_ylabel('Energy [keV]')
        self.cx.set_yscale('log')
        return

    def _plot_sst_p(self):
        """ 
        Plot the THEMIS-SST protons.
        """
        self.dx = self.fig.add_subplot(self.spec[3, :], sharex=self.bx) 
        self.dx.tick_params(axis="x", labelbottom=False) 
        dx_c = self.dx.inset_axes([1.01, 0, 0.02, 1], transform=self.dx.transAxes)
        s = SST(self.sc_id, self.time_range[0], species='i')
        s.load()
        _, p = s.spectrum(ax=self.dx, 
            pcolormesh_kwargs={'norm':matplotlib.colors.LogNorm(vmin=0.1, vmax=1E6)}
            )
        plt.colorbar(p, cax=dx_c, label='ion flux')
        self.dx.set_ylabel('Energy [keV]')
        self.dx.set_yscale('log')
        return

    def _plot_fbk(self):
        """
        Plot the THEMIS-FBK wave data
        """
        self.ex = self.fig.add_subplot(self.spec[4, :], sharex=self.bx) 
        # self.ex.tick_params(axis="x", labelbottom=False) 
        ex_c = self.ex.inset_axes([1.01, 0, 0.02, 1], transform=self.ex.transAxes)
        s = FBK(self.sc_id, self.time_range[0])
        try:
            s.load()
        except ValueError as err:
            if 'No records found for variable' in str(err):
                return
            else:
                raise
        _, p = s.spectrum(ax=self.ex, 
            pcolormesh_kwargs={'norm':matplotlib.colors.LogNorm(vmin=1E-4, vmax=0.1)})
        plt.colorbar(p, cax=ex_c, label=f'$nT^{{{2}}}/Hz$')
        self.ex.set_xlabel('Time [HH:MM]')
        self.ex.set_ylabel('Frequency [Hz]') 
        self.ex.set_yscale('log')
        self.ex.set_ylim(1, 1E3)
        time_format = matplotlib.dates.DateFormatter('%H:%M')
        self.ex.xaxis.set_major_formatter(time_format)
        return


class Sampex_Themis_ASI(Themmis_Themis_ASI):
    def __init__(self, c_filename, time_window_sec=3600/2, n_images=3, map_alt=110) -> None:
        super().__init__(c_filename, time_window_sec=time_window_sec, n_images=n_images, map_alt=map_alt)
        return

    def loop(self):
        """
        Loop over every conjunction and load the data.
        """
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


if __name__ == '__main__':
    filename = 'sampex_themis_asi_themis_aurorax_conjunctions_500_km.xlsx'
    s = Themis_Themis_ASI(filename)
    s.loop()