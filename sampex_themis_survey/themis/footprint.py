"""
Calculate the THEMIS probe footprint.
"""
import dateutil.parser
import pathlib

import IRBEM
import numpy as np
import pandas as pd

import sampex
# import pyspedas

class THEMIS_footprint:
    def __init__(self, probe, time_range):
        """"
        Load SAMPEX attitude data and calculate its footprint.

        Parameters
        ----------
        day: datetime or str
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
        self.state_files = self._get_state()
        return

    def map(self, alt=110, hemi_flag=1):
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
        m = IRBEM.MagFields(kext='OPQ77')
        _all = np.zeros_like(self.attitude.loc[:, ['Altitude', 'GEO_Lat', 'GEO_Long']])

        for i, (time, row) in enumerate(self.attitude.iterrows()):
            X = {'Time':time, 'x1':row['Altitude'], 'x2':row['GEO_Lat'], 'x3':row['GEO_Long']}
            _all[i, :] = m.find_foot_point(X, {}, alt, hemi_flag)['XFOOT']
        _all[_all == -1E31] = np.nan
        self.attitude.loc[:, ['Altitude', 'GEO_Lat', 'GEO_Long']] = _all
        return self.attitude

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
            f'th{self.probe}/l{self.level}/state/'
            f'{self.time_range[0].year}/'
                    )
        downloader = sampex.Downloader(
            base_url, 
            download_dir=self.data_dir
            )
        matched_files = downloader.ls(self.state_file_pattern)
        return