# PIRATE
**Probabilistic InteRpretation of Altimeter and in-siTu obsErvations**

This is a repo with codes for my contribution to the PIRATE project (OST-ST, 2017-2021) on work package entitled:
*Which observations/locations/timescales are most affected by Low-Frequency Chaotic Intrinsic Variability ?*

## Content
``analysis`` contains example of scripts to analyse the OCCIPUT data.

``src`` contains the ``pirate`` python package to be imported by your scripts.

## Usage
Create a dataset where we aggregate all members for a small region and a given month:
```python
import os, glob
import pirate as pr

# List all member files from 1 month:
src_path = "../data/OBS.y2014m06/"
flist = glob.glob(os.path.join(src_path + "ORCA025.L75-OCCITENS.*_enact_fdbk.nc"))

# Create and collect data of the ensemble:
ens = pr.agg.ensemble(flist, lon=[-85, -50], lat=[20, 35])
ds = ens.collect() # This can take a while for a large region
```

and plot some profile ensemble statistics:

```python
i_obs = 1
pr.plot.ensemble_profiles(ds, i_obs, vnames=['POTM', 'PSAL'])
```
![Local Profiles](/analysis/PROF_example.png?raw=true "Ensemble Profiles")

```python
i_obs = 3
pr.plot.ensemble_profile_pdf(ds, i_obs, N=15, vname="POTM", obscolor='r')
```
![Local PDF](/analysis/PDF_example.png?raw=true "Ensemble PDF")

## Sample data
One sample of data can be downloaded from here:

https://storage.googleapis.com/pirate_data/ENS/y2014m06/ORCA025.L75-OCCITENS.y2014m06_EDWregion.nc

## Dependencies
The following modules are required: ```numpy, pandas, xarray```