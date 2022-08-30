import pathlib
from datetime import datetime, timedelta

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors

import sampex
import cdflib


class FBK:
    def __init__(self, probe, day):
        """"
        Load THEMIS FBK data.

        Parameters
        ----------
        probe: str
            What 
        day: datetime
            The day to load data
        """
        self.probe = probe.lower()
        self.day = day
        self.date_str = self.day.strftime('%Y%m%d')
        self.level = 2
        self.file_pattern = (
            f'th{self.probe}_l{self.level}_fbk_{self.date_str}*.cdf'
            )
        self.data_dir = pathlib.Path(pathlib.Path.home(), 'themis-data', 
            f'th{self.probe}', 'fbk')

        self.fbk_path = self._get_fbk()
        return

    def load(self):
        """
        Load the CDF file.
        """
        self.fbk = cdflib.CDF(self.fbk_path)  # FYI self.fbk.cdf_info()['zVariables']

        self.time = np.array(self.fbk.varget(f'th{self.probe}_fb_scm1_time'), 
            dtype='datetime64[s]')
        return self.time, self.fbk

    def spectrum(self, ax=None, pcolormesh_kwargs={}):
        if ax is None:
            ax = plt.subplot()

        f = self.fbk.varget(f'th{self.probe}_fb_yaxis')
        scm1 = self.fbk.varget(f'th{self.probe}_fb_scm1')
        p = ax.pcolormesh(self.time, f, scm1.T, **pcolormesh_kwargs)
        return ax, p
        
    def _get_fbk(self):
        """
        Looks for a local fbk cdf file and downloads it if not not found.
        """
        paths = list(self.data_dir.rglob(self.file_pattern))
        if len(paths) == 0:
            # Look online
            path = self._download_fbk()
        else:
            path = paths[0]
        return path

    def _download_fbk(self):
        """
        Downloads the THEMIS probe data.

        Beware: This will only work if self.time_range spans one day.
        """
        base_url = (
            f'http://themis.ssl.berkeley.edu/data/themis/'
            f'th{self.probe}/l{self.level}/fbk/'
            f'{self.day.year}/'
                    )
        downloader = sampex.Downloader(
            base_url, 
            download_dir=self.data_dir
            )
        matched_files = downloader.ls(self.file_pattern)
        # Find the newest version
        sorted_files = sorted(matched_files, key=lambda d: d.name())
        return sorted_files[-1].download()

if __name__ == '__main__':
    s = FBK('C', datetime(2008, 3, 4, 4, 0, 0))
    s.load()
    s.spectrum(pcolormesh_kwargs={'norm':matplotlib.colors.LogNorm(vmin=1E-4, vmax=0.1)})
    plt.show()