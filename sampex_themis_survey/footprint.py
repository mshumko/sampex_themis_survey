"""
Calculate the SAMPEX footprint.
"""
import dateutil.parser

import IRBEM
import numpy as np
import pandas as pd

import sampex

class SAMPEX_footprint:
    def __init__(self, day):
        """"
        Load SAMPEX attitude data and calculate its footprint.

        Parameters
        ----------
        day: datetime or str
        """
        if isinstance(day, str):
            day = dateutil.parser.parse(day)
        day = pd.Timestamp(day).replace(hour=0, minute=0, second=0, microsecond=0)
        self.attitude = sampex.Attitude(day).load()

        self.attitude = self.attitude.loc[
            (self.attitude.index >= day) & (self.attitude.index < day+pd.Timedelta(days=1)),
            :]
        return

    def map(self, alt=110, hemi_flag=0):
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