# -*coding: UTF-8 -*-
#
# Load a dataset where we aggregated all members with "create_ensemble_for_subregion.py"
# and plot superimposed and pdf for all observations
#
__author__ = 'gmaze@ifremer.fr'

import os

wdir = os.path.expanduser('~/work/Projects/PIRATE')
# sys.path.append(os.path.join(wdir, 'src'))
os.chdir(os.path.join(wdir, 'analysis'))
import numpy as np
import xarray as xr
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt

def read_ncjuld(d):
    """Convert the netcdf Julian date into a proper format"""
    import jdcal
    from dateutil import parser
    t_ref = parser.parse(np.unique(d['JULD_REFERENCE'])[0])
    t = jdcal.jd2gcal(sum(jdcal.gcal2jd(t_ref.year, t_ref.month, t_ref.day)), np.unique(d['JULD'].values))
    return t

def ensemble_profiles(dsub, i_obs, vname = "POTM"):
    """Plot superimposed ensemble profiles
    """
    # WMO of the platform:
    wmo = np.asmatrix(dsub['STATION_IDENTIFIER'].isel(N_OBS=i_obs).values)[0, 0].strip()
    # Datetime of the measurement:
    jd = np.asmatrix(dsub['JULD'].isel(N_OBS=i_obs).values)[0, 0]
    ts = pd.Timestamp('1950-01-01') + pd.Timedelta(days=jd)
    # ts.to_datetime().strftime("%Y%m%d%H%M%S")

    vname_ens = ("%s_Hx") % (vname)
    vname_obs = ("%s_OBS") % (vname)

    y = dsub['DEPTH'].isel(N_OBS=i_obs)
    #     x_obs = dsub['PSAL_OBS'].isel(member=0).isel(N_OBS=i_obs)
    x_obs = dsub[vname_obs].isel(N_MEMBER=0).isel(N_OBS=i_obs)
    #     x_mod = dsub['PSAL_Hx'].isel(N_OBS=i_obs)
    x_mod = dsub[vname_ens].isel(N_OBS=i_obs)

    # plt.figure(dpi=100)
    sns.set(style="whitegrid")
    fig, ax = plt.subplots(nrows=1, ncols=1, figsize=(5,7))

    N = len(dsub.N_MEMBER)
    for member in np.arange(N):
        #     print dsub['PSAL_LEVEL_QC'].sel(member=member).sel(N_OBS=i_obs).values
        #         plt.plot(x_mod.sel(member=member), y, 'r', linewidth=1)
        ax.plot(x_mod.sel(N_MEMBER=member), y, '0.8', linewidth=1)
    hm, = ax.plot(x_mod.mean('N_MEMBER'), y, 'k', label='Model Ensemble mean')
    ax.plot(x_mod.mean('N_MEMBER') - x_mod.std('N_MEMBER'), y, 'k--')
    ax.plot(x_mod.mean('N_MEMBER') + x_mod.std('N_MEMBER'), y, 'k--')
    ho, = ax.plot(x_obs, y, 'r', label='Observation')
    # ax.legend([hm, ho], loc=4)
    ax.invert_yaxis()
    ax.set_xlabel(("%s [%s]")%(dsub[vname_obs].attrs['long_name'],dsub[vname_obs].attrs['units']))
    ax.set_ylabel('Depth (m)')
    # plt.title(("N_OBS #%i") % (i_obs))
    plt.title(("Measurement date: %s\nWMO: %s") % (ts.to_datetime().strftime("%Y/%m/%d %H:%M"), wmo))
    return fig

def ensemble_profile_pdf(dsub, i_obs, N=5, vname="POTM", dpi=80, hspace=-.35):
    """Plot the ensemble pdf of a vertical profile
    N is the number of depth levels
    vname = "POTM"
    vname = "PSAL"
    """
    # WMO of the platform:
    wmo = np.asmatrix(dsub['STATION_IDENTIFIER'].isel(N_OBS=i_obs).values)[0, 0].strip()
    # Datetime of the measurement:
    jd = np.asmatrix(dsub['JULD'].isel(N_OBS=i_obs).values)[0, 0]
    ts = pd.Timestamp('1950-01-01') + pd.Timedelta(days=jd)
    # ts.to_datetime().strftime("%Y%m%d%H%M%S")


    d = dsub.isel(N_OBS=i_obs)
    d = d.where(~np.isnan(d['DEPTH']), drop=True)

    iz = np.unique(np.round(np.linspace(d['N_LEVELS'].min().values,
                                        d['N_LEVELS'].max().values, N))).astype('int')

    d = d.isel(N_LEVELS=iz).stack(row=('N_LEVELS', 'N_MEMBER'))
    d = d.where(~np.isnan(d.PSAL_Hx), drop=True)
    d = d.where(~np.isnan(d.POTM_Hx), drop=True)
    df = pd.DataFrame(data={'member': d['N_MEMBER'],
                            'DEPTH': d['DEPTH'],
                            'PSAL_OBS': d['PSAL_OBS'],
                            'PSAL_Hx': d['PSAL_Hx'],
                            'POTM_OBS': d['POTM_OBS'],
                            'POTM_Hx': d['POTM_Hx']})

    # print list(df)
    # print df.index
    # print df.sample(3)

    print np.unique(d['STATION_IDENTIFIER'])
    print np.unique(d['STATION_TYPE'])
    print np.unique(d['LATITUDE'])
    print np.unique(d['LONGITUDE'])
    print np.unique(d['JULD'])
    print np.unique(d['JULD_REFERENCE'])[0]
    print read_ncjuld(d)

    print "Plot the PDF"

    vname_pdf = ("%s_Hx") % (vname)
    vname_obs = ("%s_OBS") % (vname)

    # plt.figure(dpi=300)
    sns.set(style="white", rc={"axes.facecolor": (0, 0, 0, 0)})
    # sns.set(style="dark", rc={"axes.facecolor": (0, 0, 0, 0)})

    # Initialize the FacetGrid object
    pal = sns.cubehelix_palette(N, rot=-.05, light=.7)
    # g = sns.FacetGrid(df, row="DEPTH", hue="DEPTH", aspect=5, size=1, sharey=False)
    g = sns.FacetGrid(df, row="DEPTH", hue="DEPTH", aspect=7, size=1, palette=pal, sharey=False)
    # g = sns.FacetGrid(df, row="DEPTH", hue="DEPTH", aspect=5, size=1, palette="Reds", sharey=False)
    # g = sns.FacetGrid(df, row="DEPTH", hue="DEPTH", aspect=5, size=1, palette="RdBu", sharey=False)

    g.fig.set_dpi(dpi)

    # Draw the densities in a few steps
    # g.map(sns.kdeplot, vname_pdf, clip_on=False, shade=True, alpha=1, lw=1.5, bw=.2)
    # g.map(sns.kdeplot, vname_pdf, clip_on=False, color="w", lw=2, bw=.2)
    # g.map(sns.kdeplot, vname_pdf, cumulative=True, clip_on=False, color="k", lw=2, bw=.2)
    # g.map(sns.distplot, vname_pdf, kde=True, color="b")
    g.map(sns.distplot, vname_pdf, hist=True, color="k")

    # Define and use a simple function to label the plot in axes coordinates
    def label(x, color, label):
        ax = plt.gca()
        ax.text(0, .2, ("z=%0.1fm") % (float(label)), fontweight="bold", color=color,
                ha="left", va="center", transform=ax.transAxes)
    g.map(label, vname_pdf)

    # Add observations
    def obsline(x, color, label):
        ax = plt.gca()
        plt.axvline(np.unique(x), np.min(ax.get_ylim()), np.max(ax.get_ylim()), color=color)
    g.map(obsline, vname_obs)

    # Add the horizontal 0-line:
    g.map(plt.axhline, y=0, lw=2, clip_on=False)

    # Fix axis labels:
    g.set_axis_labels(("%s [%s]") % (d[vname_pdf].attrs['long_name'], d[vname_pdf].attrs['units']))

    # Set the subplots to overlap
    g.fig.subplots_adjust(hspace=hspace)

    # Remove axes details that don't play will with overlap
    g.set_titles("")
    g.set(yticks=[])
    g.despine(bottom=True, left=True)

    # Finish
    # g.fig.suptitle("Ensemble PDF")
    g.fig.suptitle(("%s\nMeasurement: %s\nWMO: %s") % (dsub['POTM_OBS'].attrs['long_name'],
                                 ts.to_datetime().strftime("%Y/%m/%d %H:%M"), wmo))

    return g


print "Load the dataset where we aggregated all members for a small region"
# dsub = xr.open_dataset('tempo.nc', decode_times=False, mask_and_scale=True, cache=True)
dsub = xr.open_dataset('ORCA025.L75-OCCITENS-smallEDWregion.nc', decode_times=False, mask_and_scale=True, cache=True)

# Super imposed model data (gray) and observation profile (red), and PDF:
print "Plot all profiles:"
vname = "POTM"

obs = np.arange(len(dsub.N_OBS))
# obs = [2,4,11,12,15]
# obs = [2]
for i_obs in obs:
    # Misc info about the obs:
    wmo = np.asmatrix(dsub['STATION_IDENTIFIER'].isel(N_OBS=i_obs).values)[0, 0].strip()
    ts = pd.Timestamp('1950-01-01') \
         + pd.Timedelta(days=np.asmatrix(dsub['JULD'].isel(N_OBS=i_obs).values)[0, 0])
    date_str = ts.to_datetime().strftime("%Y%m%d%H%M")

    # Plot 1
    f = ensemble_profiles(dsub, i_obs, vname = vname)
    # Save it:
    f.savefig(os.path.join("img", ("%s_%s_%s.pdf")%(date_str, wmo, vname)))

    # Plot 2:
    g = ensemble_profile_pdf(dsub, i_obs, N=15, vname=vname)
    # Save it:
    g.savefig(os.path.join("img", ("%s_%s_%s_pdf.pdf")%(date_str, wmo, vname)))

    #
    plt.close('all')