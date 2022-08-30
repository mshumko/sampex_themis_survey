"""
Calculate the THEMIS probe footprint.
"""
import dateutil.parser
import pathlib
from datetime import datetime

import IRBEM
import numpy as np
import pandas as pd

import sampex
# import pyspedas
import cdflib

re_km = 6371

class THEMIS_footprint:
    def __init__(self, probe, time_range):
        """"
        Load THEMIS state data and calculate its footprint.

        Parameters
        ----------
        probe: str
            What 
        time_range: datetime
            The data time range
        """
        self.probe = probe
        self.time_range = time_range
        self.date_str = self.time_range[0].strftime('%Y%m%d')
        self.level = 1
        self.state_file_pattern = (
            f'th{self.probe.lower()}_l{self.level}_state_{self.date_str}*.cdf'
            )
        self.data_dir = pathlib.Path(pathlib.Path.home(), 'themis-data', 
            f'th{self.probe.lower()}', 'state')

        self.state_path = self._get_state()
        self.state = self._load_state()
        return

    def map(self, alt=110, hemi_flag=1, model='OPQ77', mag_input={}):
        """
        Map self.lla along the magnetic field line to alt using IRBEM.MagFields.find_foot_print.

        Parameters
        ----------
        alt: float
            The mapping altitude in units of kilometers
        hemi_flag: int
            What direction to trace the field line: 
            0 = same magnetic hemisphere as starting point
            +1   = northern magnetic hemisphere
            -1   = southern magnetic hemisphere
            +2   = opposite magnetic hemisphere as starting point
        """
        m = IRBEM.MagFields(kext=model)
        _all = np.zeros_like(self.lla)

        for i, (time, coord) in enumerate(zip(self.time, self.lla)):
            X = {'Time':time, 'x1':coord[0], 'x2':coord[1], 'x3':coord[2]}
            _all[i, :] = m.find_foot_point(X, mag_input, alt, hemi_flag)['XFOOT']
        _all[_all == -1E31] = np.nan
        self.footprint = np.roll(_all, -1, axis=1)
        return self.footprint

    def _get_state(self):
        """
        Looks for a local state file and downloads it if not not found.
        """
        paths = list(self.data_dir.rglob(self.state_file_pattern))
        if len(paths) == 0:
            # Look online
            path = self._download_state()
        else:
            path = paths[0]
        return path

    def _download_state(self):
        """
        Downloads the THEMIS probe data.

        Beware: This will only work if self.time_range spans one day.
        """
        base_url = (
            f'http://themis.ssl.berkeley.edu/data/themis/'
            f'th{self.probe.lower()}/l{self.level}/state/'
            f'{self.time_range[0].year}/'
                    )
        downloader = sampex.Downloader(
            base_url, 
            download_dir=self.data_dir
            )
        matched_files = downloader.ls(self.state_file_pattern)
        # Find the newest version
        sorted_files = sorted(matched_files, key=lambda d: d.name())
        return sorted_files[-1].download()

    def _load_state(self):
        """
        Load the CDF file.
        """
        self.state = cdflib.CDF(self.state_path)

        self.time = np.array(self.state.varget(f'th{self.probe.lower()}_state_time'), 
            dtype='datetime64[s]')
        # Convert to datetime.datetime for IRBEM
        self.time = pd.to_datetime(self.time).to_pydatetime()
        self.gei = self.state.varget(f'th{self.probe.lower()}_pos')/re_km

        # Time filter
        valid_idt = np.where(
            (self.time >= self.time_range[0]) & 
            (self.time <= self.time_range[1])
            )[0]
        self.time = self.time[valid_idt]
        self.gei = self.gei[valid_idt, :]
        
        self.lla = IRBEM.Coords().transform(self.time, self.gei, 'GEI', 'GDZ')
        return 