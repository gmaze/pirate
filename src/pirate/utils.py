# -*coding: UTF-8 -*-
__author__ = 'gmaze@ifremer.fr'

import jdcal
from dateutil import parser
import numpy as np

def read_ncjuld(d):
    """Convert the netcdf Julian date into a proper format"""
    t_ref = parser.parse(np.unique(d['JULD_REFERENCE'])[0])
    t = jdcal.jd2gcal(sum(jdcal.gcal2jd(t_ref.year, t_ref.month, t_ref.day)),np.unique(d['JULD']))
    return t

# def open_erddap(dsub, i_obs):
#     # WMO of the platform:
#     wmo = np.asmatrix(dsub['STATION_IDENTIFIER'].isel(N_OBS=i_obs).values)[0, 0].strip()
#     # Datetime of the measurement:
#     jd = np.asmatrix(dsub['JULD'].isel(N_OBS=i_obs).values)[0, 0]
#     ts = pd.Timestamp('1950-01-01') + pd.Timedelta(days=jd)
#     # ts.to_datetime().strftime("%Y%m%d%H%M%S")
#     # http://www.ifremer.fr/erddap/tabledap/ArgoFloats.htmlTable?temp%2Cpres%2Cpsal%2Ctime&platform_number=%226901047%22&time=2014-06-07T05%3A53%3A50Z