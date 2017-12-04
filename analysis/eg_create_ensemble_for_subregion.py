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
import numpy as np
import xarray as xr

# Define where to get the data for a single month:
# src_path = "../data/sample/OBS.y1993m01/"
src_path = "../data/OCCOBS/OBS.y2014m06/"
flist = glob.glob(os.path.join(src_path + "ORCA025.L75-OCCITENS.*_enact_fdbk.nc"))
print flist

# Define functions for pre-processing of each member file:

def subselect_region(ds):
    #     mask = dsub['PSAL_OBS'] != dsub['PSAL_OBS']._Fillvalue
    #     ds.coords['mask'] = (('member', 'N_OBS', 'N_LEVELS'), mask)
    #     ds = ds.where(ds.mask >= 1, drop=True)
    ds = ds.where((ds.LATITUDE >= 20) & (ds.LATITUDE <= 45) & (ds.LONGITUDE >= -85) & (ds.LONGITUDE <= -30), drop=True)
    return ds


def subselect_qc(dsub):
    mask = (np.abs(dsub['POTM_OBS']) != dsub['POTM_OBS']._Fillvalue) \
            & (dsub['POTM_Hx'] != dsub['POTM_Hx']._Fillvalue) \
            & (np.abs(dsub['PSAL_OBS']) != dsub['PSAL_OBS']._Fillvalue) \
            & (dsub['PSAL_Hx'] != dsub['PSAL_Hx']._Fillvalue) \
            & (dsub['POSITION_QC'] == 1.) \
            & (dsub['DEPTH_QC'] == 1.) \
            & (dsub['PSAL_LEVEL_QC'] == 1.) \
            & (dsub['POTM_LEVEL_QC'] == 1.)
    dsub.coords['mask'] = (('N_OBS', 'N_LEVELS'), mask)
    # dsub.coords['LONGITUDE'] = (('N_OBS'), dsub.LONGITUDE)
    # dsub.coords['LATITUDE'] = (('N_OBS'), dsub.LATITUDE)
    # dsub.coords['DEPTH'] = (('N_OBS', 'N_LEVELS'), dsub.DEPTH)
    # dsub.coords['JULD'] = (('N_OBS'), dsub.JULD)
    # dsub.coords['STATION_IDENTIFIER'] = (('N_OBS'), dsub.STATION_IDENTIFIER)
    # dsub.coords['STATION_TYPE'] = (('N_OBS'), dsub.STATION_TYPE)
    dsub = dsub.where(dsub.mask >= 1, drop=True)
    dsub = dsub.drop('mask')
    return dsub


def subselect_Nmin(ds):
    """"Remove profiles with not enough data
        We keep only profiles with at least 10 values over a 1000m thickness layer
    """
    # Nb of valid data for each profile:
    # len(ds['N_LEVELS'])-ds['POTM_OBS'].isnull().sum('N_LEVELS')
    N = ds['POTM_OBS'].notnull().sum('N_LEVELS')
    ds['N'] = N

    # Thickness of the valid points layer:
    H = ds['DEPTH'].where(ds['POTM_OBS'].notnull()).max(dim='N_LEVELS')-\
    ds['DEPTH'].where(ds['POTM_OBS'].notnull()).min(dim='N_LEVELS')
    ds['H'] = H

    ds['KEEP'] = xr.DataArray(
        np.all((ds['H'] >= 1000, ds['N'] >= 10), axis=0),
        dims= {'N_OBS':ds['N_OBS']})

    ds = ds.where(ds['KEEP'], drop=True)
    ds = ds.drop(['N', 'H', 'KEEP'])

    return ds


def subselect(ds):
    """ Run all selection methods"""
    ds = subselect_region(ds)
    ds = subselect_qc(ds)
    ds = subselect_Nmin(ds)
    return ds


# Create the ensemble dataset:
ds = xr.open_mfdataset(flist, decode_times=False, concat_dim='N_MEMBER', preprocess=subselect, mask_and_scale=True)
# ds.coords['LATITUDE'] = (('N_MEMBER', 'N_OBS', 'N_LEVELS'), ds.LATITUDE)
# ds.coords['LONGITUDE'] = (('N_MEMBER', 'N_OBS', 'N_LEVELS'), ds.LONGITUDE)
# ds.coords['DEPTH'] = (('N_MEMBER', 'N_OBS', 'N_LEVELS'), ds.DEPTH)
# ds.coords['JULD'] = (('N_MEMBER', 'N_OBS', 'N_LEVELS'), ds.JULD)
# ds.coords['STATION_IDENTIFIER'] = (('N_MEMBER', 'N_OBS', 'N_LEVELS'), ds.STATION_IDENTIFIER)
# ds.coords['STATION_TYPE'] = (('N_MEMBER', 'N_OBS', 'N_LEVELS'), ds.STATION_TYPE)

# Now remove unnecessary dimensions to some variables
# Redundant information through the N_MEMBER dimension:
vlist = ['VARIABLES', 'ENTRIES', 'EXTRA', 'DEPTH_QC', 'DEPTH_QC_FLAGS',
         'JULD_REFERENCE', 'OBSERVATION_QC', 'OBSERVATION_QC_FLAGS',
         'POSITION_QC', 'POSITION_QC_FLAGS', 'JULD_QC', 'JULD_QC_FLAGS',
         'ORIGINAL_FILE_INDEX', 'LATITUDE', 'LONGITUDE', 'JULD',
         'STATION_IDENTIFIER', 'STATION_TYPE', 'DEPTH']
for v in vlist:
    ds[v] = ds[v].isel(N_MEMBER=0)

# Redundant information through the N_LEVELS dimension:
vlist = ['LATITUDE', 'LONGITUDE', 'JULD', 'VARIABLES',
         'STATION_IDENTIFIER', 'STATION_TYPE',
         'ENTRIES', 'EXTRA', 'JULD_REFERENCE',
         'POSITION_QC', 'POSITION_QC_FLAGS', 'JULD_QC_FLAGS']
for v in vlist:
    ds[v] = ds[v].isel(N_LEVELS=0)

# Set some variables as coordinates:
ds.coords['LATITUDE'] = (('N_OBS'), ds.LATITUDE)
ds.coords['LONGITUDE'] = (('N_OBS'), ds.LONGITUDE)
ds.coords['JULD'] = (('N_OBS'), ds.JULD)
ds.coords['DEPTH'] = (('N_OBS', 'N_LEVELS'), ds.DEPTH)

print ds

# Load data into memory:
# ds.compute() # This is very long !

# Save it on disk (shorter if compute was executed before):
# ds.to_netcdf('ORCA025.L75-OCCITENS-smallEDWregion.nc')  # This is very long !
ds.to_netcdf('ORCA025.L75-OCCITENS-EDWregion.nc')  # This is very long !