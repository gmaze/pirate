# -*coding: UTF-8 -*-
__author__ = 'gmaze@ifremer.fr'

import getpass

import datetime
import numpy as np
import xarray as xr


class ensemble:
    """
        import pirate as pr
        flist = glob("*")
        ens = pr.agg.ensemble(flist, lon=[-85, -30], lat=[20, 45])
        ds = ens.collect()
    """
    def __init__(self, flist, lon=[-180, 180], lat=[-90, 90], Nmin=10, Hmin=1000):
        self._props = {'flist': flist,
                       'lon': lon,
                       'lat': lat,
                       'Nmin': Nmin,
                       'Hmin': Hmin}

    def collect(self):
        """Return the dataset with the ensemble sub-selected"""

        # Create the subselecter functions:
        def subselect(ds, props):
            """ chain all selection methods
                ds = subselect(ds, lon=[-180, 180], lat=[-90, 90], Nmin=10, Hmin=1000)
            """
            def subselect_Nmin(ds, props):
                """"Remove profiles with not enough data
                    We keep only profiles with at least 10 values over a 1000m thickness layer
                    ds = subselect_Nmin(ds, Nmin=10, Hmin=1000)
                """
                Nmin = props['Nmin']
                Hmin = props['Hmin']

                # Nb of valid data for each profile:
                # len(ds['N_LEVELS'])-ds['POTM_OBS'].isnull().sum('N_LEVELS')
                N = ds['POTM_OBS'].notnull().sum('N_LEVELS')
                # print N
                ds['N'] = N

                # Thickness of the valid points layer:
                H = ds['DEPTH'].where(ds['POTM_OBS'].notnull()).max(dim='N_LEVELS') - \
                    ds['DEPTH'].where(ds['POTM_OBS'].notnull()).min(dim='N_LEVELS')
                # print H.values
                ds['H'] = H

                ds['KEEP'] = xr.DataArray(
                    np.all((ds['H'] >= Hmin, ds['N'] >= Nmin), axis=0),
                    dims={'N_OBS': ds['N_OBS']})

                ds = ds.where(ds['KEEP'], drop=True)
                ds = ds.drop(['N', 'H', 'KEEP'])

                return ds

            def subselect_qc(dsub):
                """Drop profiles and measurements with bad Quality Flags
                    ds = subselect_qc(ds)
                """
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

            def subselect_region(ds, props):
                """Select profiles for a given rage of latitude/longitude
                    ds = subselect_region(ds, lon=[-85, -30], lat=[20, 45])
                """
                # ds = ds.where((ds.LATITUDE >= 20) & (ds.LATITUDE <= 45) & (ds.LONGITUDE >= -85) & (ds.LONGITUDE <= -30), drop=True)
                ds = ds.where(
                    (ds.LATITUDE >= props['lat'][0]) & (ds.LATITUDE <= props['lat'][1]) &
                    (ds.LONGITUDE >= props['lon'][0]) & (ds.LONGITUDE <= props['lon'][1]),
                    drop=True)
                return ds

            ds = subselect_region(ds, props)
            ds = subselect_qc(ds)
            ds = subselect_Nmin(ds, props)
            return ds

        def subselecter(ds):
            props = self._props
            return subselect(ds, props)

        # Open the ensemble dataset and appy sub-selecter:
        ds = xr.open_mfdataset(self._props['flist'],
                               decode_times=False,
                               concat_dim='N_MEMBER',
                               preprocess=subselecter,
                               mask_and_scale=True)
        # print ds
        # print "- " * 24

        if len(ds['N_OBS']) == 0:
            raise Exception("No observations match your selection parameters !")

        if len(ds['N_LEVELS']) == 0:
            raise Exception("No observations match your selection parameters !")

        # Now remove unnecessary dimensions to some variables
        # Redundant information through the N_MEMBER dimension:
        vlist = ['VARIABLES', 'ENTRIES', 'EXTRA', 'DEPTH_QC', 'DEPTH_QC_FLAGS',
                 'JULD_REFERENCE', 'OBSERVATION_QC', 'OBSERVATION_QC_FLAGS',
                 'POSITION_QC', 'POSITION_QC_FLAGS', 'JULD_QC', 'JULD_QC_FLAGS',
                 'ORIGINAL_FILE_INDEX', 'LATITUDE', 'LONGITUDE', 'JULD',
                 'STATION_IDENTIFIER', 'STATION_TYPE', 'DEPTH']
        for v in vlist:
            # print "Remove N_MEMBER dimension to:", v
            ds[v] = ds[v].isel(N_MEMBER=0)

        # return ds
        # print ds
        # print ds['N_LEVELS']
        # print ds['N_OBS']
        # print ds['N_OBS'] == 0

        # print "- " * 24
        # Redundant information through the N_LEVELS dimension:
        vlist = ['LATITUDE', 'LONGITUDE', 'JULD', 'VARIABLES',
                 'STATION_IDENTIFIER', 'STATION_TYPE',
                 'ENTRIES', 'EXTRA', 'JULD_REFERENCE', 'JULD_QC',
                 'POSITION_QC', 'POSITION_QC_FLAGS', 'JULD_QC_FLAGS']
        for v in vlist:
            # print "Remove N_LEVELS dimension to:", v
            ds[v] = ds[v].isel(N_LEVELS=0)

        # print ds
        # Set some variables as coordinates:
        ds.coords['LATITUDE'] = (('N_OBS'), ds.LATITUDE)
        ds.coords['LONGITUDE'] = (('N_OBS'), ds.LONGITUDE)
        ds.coords['JULD'] = (('N_OBS'), ds.JULD)
        ds.coords['DEPTH'] = (('N_OBS', 'N_LEVELS'), ds.DEPTH)

        # Add more global attributes:
        ds.attrs['Creation date'] = datetime.datetime.now().strftime('%Y/%m/%d %H:%M:%S')
        ds.attrs['Created by'] = getpass.getuser()
        ds.attrs['Software'] = 'http://github.com/gmaze/pirate'
        ds.attrs['Domain'] = ("%0.2f<LON<%02.f and %0.2f<LAT<%02.f")%\
                                      (self._props['lon'][0], self._props['lon'][1],
                                       self._props['lat'][0], self._props['lat'][1])
        ds.attrs['Vertical restriction'] = ("Valid observations with more than %i points "
                                            "over a %0.0fm layer")%(
                                        self._props['Nmin'], self._props['Hmin'])

        return ds