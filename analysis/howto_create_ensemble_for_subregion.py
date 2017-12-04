#!/usr/bin/env python
# -*coding: UTF-8 -*-
#
# Here we create a single netcdf file for a given month of ensemble data.
# Observations are selected for a given region and with a minimum number of samples over depth.
# A netcdf file is produced with aggregated model data along a new "N_MEMBER" dimension
#
# This methodology is not efficient, it takes a very long time (hours) to
# aggregate data for a single month over a box of 25x55 deg (N_OBS = 458)
#
__author__ = 'gmaze@ifremer.fr'

import os
import sys
wdir = os.path.expanduser('~/work/Projects/PIRATE')
sys.path.append(os.path.join(wdir, 'src'))
os.chdir(os.path.join(wdir, 'analysis'))
import glob
import pirate as pr

# List member files for a single month:
src_path = "../data/sample/OBS.y1993m01/"
# src_path = "../data/OCCOBS/OBS.y2014m06/"
flist = glob.glob(os.path.join(src_path + "ORCA025.L75-OCCITENS.*_enact_fdbk.nc"))
print flist

# Create and collect the ensemble dataset data

ens = pr.agg.ensemble(flist, lon=[-85, -50], lat=[20, 35])
ds = ens.collect()
pr.plot.ensemble_profiles(ds, 0, vnames='POTM')
pr.plot.ensemble_profiles(ds, 0, vnames=['POTM', 'PSAL'])

ens = pr.agg.ensemble(flist, lon=[-55, -30], lat=[20, 35], Nmin=10)
ds = ens.collect()

ens = pr.agg.ensemble(flist, lon=[-85, -30], lat=[20, 45], Nmin=10, Hmin=1000)
ds = ens.collect()

print ds

# Load data into memory:
# ds.compute() # This is very long !

# Save it on disk (shorter if compute was executed before):
# ds.to_netcdf('ORCA025.L75-OCCITENS-smallEDWregion.nc')  # This is very long !
ds.to_netcdf('ORCA025.L75-OCCITENS-EDWregion.nc')  # This is very long !