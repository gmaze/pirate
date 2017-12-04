# -*coding: UTF-8 -*-

__author__ = 'gmaze@ifremer.fr'
import os
import sys

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import xarray as xr

wdir = os.path.expanduser('~/work/Projects/PIRATE')
os.chdir(os.path.join(wdir, 'analysis'))

# Now import our PIRATE machinery:
sys.path.append(os.path.join(wdir, 'src'))
import pirate as pr

def footprint(f):
    plt.figtext(0.01, 0.01, ("%s / PIRATE Project") % (__author__),
                figure=f, horizontalalignment='left', fontsize=6,
                rotation='vertical', rotation_mode='anchor')

#
print "Load the dataset where we aggregated all members for a small region"
# dsub = xr.open_dataset('tempo.nc', decode_times=False, mask_and_scale=True, cache=True)
dsub = xr.open_dataset('ORCA025.L75-OCCITENS-smallEDWregion.nc', decode_times=False, mask_and_scale=True, cache=True)

# Superimposed model data (gray) and observation profile (red):
print "Plot all profiles:"

obs = np.arange(len(dsub.N_OBS))
# obs = [2,4,11,12,15]
# obs = [2]
for i_obs in obs:
    # Misc info about the obs:
    wmo = np.asmatrix(dsub['STATION_IDENTIFIER'].isel(N_OBS=i_obs).values)[0, 0].strip()
    ts = pd.Timestamp('1950-01-01') \
         + pd.Timedelta(days=np.asmatrix(dsub['JULD'].isel(N_OBS=i_obs).values)[0, 0])
    date_str = ts.to_pydatetime().strftime("%Y%m%d%H%M")

    # Plots:
    vlist = ["POTM","PSAL"]
    for vname in vlist:
        f = pr.plot.ensemble_profiles(dsub, i_obs, vname = vname)
        footprint(f)
        # Save it:
        f.savefig(os.path.join("img", vname, "PROF", ("%s_%s_%s.pdf")%(date_str, wmo, vname)))

    #
    plt.close('all')


# Model PDF at selected depth levelsL
print "Plot all profiles PDF:"
vname = "POTM"

obs = np.arange(len(dsub.N_OBS))
# obs = [2,4,11,12,15]
# obs = [2]
for i_obs in obs:
    # Misc info about the obs:
    wmo = np.asmatrix(dsub['STATION_IDENTIFIER'].isel(N_OBS=i_obs).values)[0, 0].strip()
    ts = pd.Timestamp('1950-01-01') \
         + pd.Timedelta(days=np.asmatrix(dsub['JULD'].isel(N_OBS=i_obs).values)[0, 0])
    date_str = ts.to_pydatetime().strftime("%Y%m%d%H%M")

    # Plot:
    vlist = ["POTM","PSAL"]
    for vname in vlist:
        g = pr.plot.ensemble_profile_pdf(dsub, i_obs, N=15, vname=vname, obscolor='r')
        # g = pr.plot.ensemble_profile_pdf(dsub, i_obs, N=5, vname=vname, obscolor='r')
        # g = pr.plot.ensemble_profile_pdf(dsub, i_obs, N=5, hspace=0, vname=vname, obscolor='r')
        footprint(g.fig)
        # Save it:
        g.savefig(os.path.join("img", vname, "PDF", ("%s_%s_%s_pdf.pdf")%(date_str, wmo, vname)))

    #
    plt.close('all')
