"""
Download and load THEMIS-SST data
"""
import pathlib
from datetime import datetime

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors

import sampex
import cdflib


class SST:
    def __init__(self, probe, day, species='e'):
        """"
        Load THEMIS SST data.

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
        self.species = species.lower()
        assert self.species.lower() in ['e', 'i']
        self.file_pattern = (
            f'th{self.probe}_l{self.level}_sst_{self.date_str}*.cdf'
            )
        self.data_dir = pathlib.Path(pathlib.Path.home(), 'themis-data', 
            f'th{self.probe}', 'sst')

        self.sst_path = self._get_sst()
        return

    def load(self):
        """
        Load the CDF file.
        """
        self.sst = cdflib.CDF(self.sst_path)  # FYI self.sst.cdf_info()['zVariables']

        self.time = np.array(self.sst.varget(f'th{self.probe}_ps{self.species}f_time'), 
            dtype='datetime64[s]')
        return self.time, self.sst

    def spectrum(self, ax=None, pcolormesh_kwargs={}):
        if ax is None:
            ax = plt.subplot()
        E = self.sst.varget(f'th{self.probe}_ps{self.species}f_en_eflux_yaxis')[0, :]/1E3
        # Finite energies
        ide = np.where(np.isfinite(E))[0]
        E = E[ide]
        eflux = self.sst.varget(f'th{self.probe}_ps{self.species.lower()}f_en_eflux')
        eflux = eflux[:, ide]
        p = ax.pcolormesh(self.time, E, eflux.T, **pcolormesh_kwargs)
        return ax, p
        
    def _get_sst(self):
        """
        Looks for a local SST cdf file and downloads it if not not found.
        """
        paths = list(self.data_dir.rglob(self.file_pattern))
        if len(paths) == 0:
            # Look online
            path = self._download_sst()
        else:
            path = paths[0]
        return path

    def _download_sst(self):
        """
        Downloads the THEMIS probe data.

        Beware: This will only work if self.time_range spans one day.
        """
        base_url = (
            f'http://themis.ssl.berkeley.edu/data/themis/'
            f'th{self.probe}/l{self.level}/sst/'
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
    s = SST('E', datetime(2008, 3, 4, 4, 0, 0),species='i')
    s.load()
    s.spectrum(pcolormesh_kwargs={'norm':matplotlib.colors.LogNorm(vmin=0.1, vmax=1E6)})
    plt.show()