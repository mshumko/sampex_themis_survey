"""
Plot the instantaneous locations of the THEMIS ASIs, and THEMIS probe footprints 
during  the 2008-03-04 triple conjunction.
"""
from datetime import datetime, date, timedelta
from importlib.resources import path
import string
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


class Summary_Plot:
    def __init__(self, c_filename, time_window_sec=3600/2, n_images=3, map_alt=110) -> None:
        self.c_filename = c_filename
        self.time_window_sec = time_window_sec
        self.n_images = n_images
        self.map_alt = map_alt
        self._load()
        return

    def _load(self):
        self.c_path = pathlib.Path(sampex_themis_survey.config['code_dir'], 'data', self.c_filename)
        self.c_df = pd.read_excel(self.c_path, skiprows=1)
        self.c_df['start'] = pd.to_datetime(self.c_df['Start Time (UTC)'])
        self.c_df['end'] = pd.to_datetime(self.c_df['End Time (UTC)'])
        return

    def loop(self):
        """
        Loop over every conjunction and load the data.
        """
        current_date = date.min
        for _, row in self.c_df.iterrows():
            print(f'Processing {row["start"]}')
            if current_date != row['start'].date:
                try:
                    self.hilt = sampex.HILT(row['start']).load()
                except FileNotFoundError as err:
                    if 'does not contain any hyper references' in str(err):
                        continue
                    raise
                current_date = row['start'].date

            self.asi_location = row['Conjunction Between'].split('and')[0].rstrip().split(' ')[1]
            self.sc_id = row['Conjunction Between'].split('and')[1].rstrip().split('-')[1]
            self.time_range = [
                row['start'] - timedelta(seconds=self.time_window_sec/2),
                row['end'] + timedelta(seconds=self.time_window_sec/2)
            ]
            # Make the plot here
            self.plot()
        
        return

    def plot(self):
        self.fig = plt.figure(figsize=(12, 9))
        self.spec = gridspec.GridSpec(
            nrows=5, 
            ncols=self.n_images, figure=self.fig, height_ratios=(2, 1, 1, 1, 1)
            )
        # -1 so that image_time[-1] == time_range[1]
        dt = (self.time_range[1]-self.time_range[0])/(self.n_images-1)
        self.image_times = [self.time_range[0] +  i*dt for i in range(self.n_images)]
        self._plot_asi_images()
        self._plot_themis_footprint()
        self._plot_keogram()
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
            0.9*np.min(skymap['SITE_MAP_LATITUDE']),
            1.1*np.max(skymap['SITE_MAP_LATITUDE'])
            )
        self.lon_bounds = (
            0.9*np.min(skymap['SITE_MAP_LONGITUDE']),
            1.1*np.max(skymap['SITE_MAP_LONGITUDE'])
            )

        for i, (image_time, subplot_letter) in enumerate(z):
            self.ax[i] = self.fig.add_subplot(self.spec[0, i])
            self.ax[i].get_xaxis().set_visible(False)
            self.ax[i].get_yaxis().set_visible(False)
            asilib.make_map(lat_bounds=self.lat_bounds, lon_bounds=self.lon_bounds,
                ax=self.ax[i], land_color='grey')

            # Update image_times with the actual image timestamp.
            self.image_times[i], _, _, _, _ = asilib.plot_map(
                'THEMIS', self.asi_location, image_time, self.map_alt, 
                ax=self.ax[i], asi_label=False, color_bounds=None, pcolormesh_kwargs={'rasterized':True}
                )
            
            self.ax[i].text(0, 0, f'({subplot_letter}) {image_time.strftime("%H:%M:%S")}', 
                transform=self.ax[i].transAxes, va='bottom', color='green', fontsize=15)
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
                color='r', s=100, marker='.')
        return

    def _plot_keogram(self):
        """
        Plots a keogram on the 2nd row of the subplot grid.
        """
        self.bx = self.fig.add_subplot(self.spec[1, :])
        self.bx.tick_params(axis="x", labelbottom=False)
        asilib.plot_keogram('THEMIS', self.asi_location, self.time_range, 
            ax=self.bx, map_alt=self.map_alt, aacgm=True)
        self.bx.set_xlim(self.time_range)
        self.bx.set_title('')
        self.bx.text(0, 1, f'({string.ascii_uppercase[self.n_images+1]}) THEMIS ASI keogram', 
            transform=self.bx.transAxes, va='top', color='white', fontsize=15)
        self.bx.set_ylabel('Magnetic lat [geg]')
        # self.bx.set_ylim(self.lat_bounds)
        return


if __name__ == '__main__':
    s = Summary_Plot('sampex_themis_asi_themis_aurorax_conjunctions_500_km.xlsx')
    s.loop()

# themis_footprints = {}
# colors = ['r', 'c', 'b', 'k']
# for themis_probe, color in zip(themis_probes, colors):
#     f = THEMIS_footprint(themis_probe, (times[0], times[-1]))
#     f.map(alt=alt, model='T89', mag_input={'Kp':1.667})
#     themis_footprints[themis_probe] = {'time':f.time, 'footprint':f.footprint}
    
#     for time, ax_i in zip(times, ax):
#         dt = np.abs([(time-ti).total_seconds() for ti in f.time])
#         idt = np.argmin(dt)
#         assert (time-f.time[idt]).total_seconds() <= 60, (
#             f'The time stamp difference, {(time-f.time[idt]).total_seconds()} '
#             f's, is too large.')
#         ax_i.scatter(f.footprint[idt, 1], f.footprint[idt, 0], 
#             label=f'THEMIS-{themis_probe.upper()}',
#             color=color, s=100, marker='.')

# for i, themis_asi_location in enumerate(themis_asi_locations, start=1):
#     bx = fig.add_subplot(spec[i, :])
#     bx.tick_params(axis="x", labelbottom=False)
#     asilib.plot_keogram('THEMIS', themis_asi_location, plot_range, 
#         ax=bx, map_alt=alt, aacgm=True)
#     bx.set_xlim(times[0], times[-1])
#     bx.set_title('')
#     bx.text(0, 1, f'({string.ascii_uppercase[n+i-1]}) THEMIS ASI keogram', 
#         transform=bx.transAxes, va='top', color='white', fontsize=15)
#     bx.set_ylabel('Magnetic lat [geg]')
#     bx.set_ylim(66, 72)

# cx = [fig.add_subplot(spec[i, :], sharex=bx) 
#     for i in range(2, len(themis_asi_locations)+4)]
# for cx_i in cx[:-1]:
#     cx_i.tick_params(axis="x", labelbottom=False) 

# # Electrons
# cx_c = cx[0].inset_axes([1.01, 0, 0.02, 1], transform=cx[0].transAxes)
# s = SST(themis_probes[0], plot_range[0], species='e')
# s.load()
# _, p = s.spectrum(ax=cx[0], 
#     pcolormesh_kwargs={'norm':matplotlib.colors.LogNorm(vmin=1E2, vmax=5E5)}
#     )
# plt.colorbar(p, ax=cx, cax=cx_c, label='Electron flux')
# cx[0].set_xlim(plot_range)
# cx[0].set_ylabel('Energy [keV]')
# cx[0].set_yscale('log')

# # Protons
# cx_c = cx[1].inset_axes([1.01, 0, 0.02, 1], transform=cx[1].transAxes)
# s = SST(themis_probes[0], plot_range[0], species='i')
# s.load()
# _, p = s.spectrum(ax=cx[1], 
#     pcolormesh_kwargs={'norm':matplotlib.colors.LogNorm(vmin=0.1, vmax=1E6)}
#     )
# plt.colorbar(p, ax=cx[1], cax=cx_c, label='ion flux')
# cx[1].set_xlim(plot_range)
# cx[1].set_ylabel('Energy [keV]')
# cx[1].set_yscale('log')

# # Waves
# s = FBK(themis_probes[0], plot_range[0])
# s.load()
# s.spectrum(ax=cx[2], 
#     pcolormesh_kwargs={'norm':matplotlib.colors.LogNorm(vmin=1E-4, vmax=0.1)})
# cx[2].set_xlabel('Time [HH:MM]')
# cx[2].set_ylabel('Frequency [Hz]') 
# cx[2].set_yscale('log')
# cx[2].set_ylim(1, 1E3)
# time_format = matplotlib.dates.DateFormatter('%H:%M')
# cx[2].xaxis.set_major_formatter(time_format)

# plot_labels = (
#     'THEMIS-SST electrons',
#     'THEMIS-SST ions',
#     'THEMIS-FBK'
# )
# themis_label_y = (0, 0, 0.8)

# for i, (cx_i, plot_label, y) in enumerate(zip(cx, plot_labels, themis_label_y), start=n+1):
#     cx_i.text(0, y, f'({string.ascii_uppercase[i]}) {plot_label}', 
#         transform=cx_i.transAxes, va='bottom', color='white', fontsize=15)

# plt.suptitle(f'Example Triple Conjunction | THEMIS-THEMIS ASI-SAMPEX', fontsize=15)

# # plt.legend()
# plt.subplots_adjust(wspace=0.02, hspace=0.07, left=0.055, right=0.92, top=0.943)
# plt.show()