# -*coding: UTF-8 -*-
__author__ = 'gmaze@ifremer.fr'

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns

import pirate as pr


def ensemble_profiles(dsub, i_obs, vnames = "POTM"):
    """Plot ensemble profiles superimposed with observation profile
        vname = "POTM"
        vname = "PSAL"
    """
    # WMO of the platform:
    wmo = np.asmatrix(dsub['STATION_IDENTIFIER'].isel(N_OBS=i_obs).values)[0, 0].strip()
    # Datetime of the measurement:
    jd = np.asmatrix(dsub['JULD'].isel(N_OBS=i_obs).values)[0, 0]
    ts = pd.Timestamp('1950-01-01') + pd.Timedelta(days=jd)
    # ts.to_pydatetime().strftime("%Y%m%d%H%M%S")

    if len(vnames) == 2:
        fig, ax = plt.subplots(nrows=1, ncols=2, figsize=(10,7), sharey=True)
    else:
        vnames = [vnames]
        fig, ax = plt.subplots(nrows=1, ncols=1, figsize=(5,7))
        ax = np.array([ax])

    # plt.figure(dpi=100)
    sns.set(style="whitegrid")
    N = len(dsub.N_MEMBER)

    for vname, iv in zip(vnames, range(len(vnames))):
        vname_ens = ("%s_Hx") % (vname)
        vname_obs = ("%s_OBS") % (vname)

        y = dsub['DEPTH'].isel(N_OBS=i_obs)
        #     x_obs = dsub['PSAL_OBS'].isel(member=0).isel(N_OBS=i_obs)
        x_obs = dsub[vname_obs].isel(N_MEMBER=0).isel(N_OBS=i_obs)
        #     x_mod = dsub['PSAL_Hx'].isel(N_OBS=i_obs)
        x_mod = dsub[vname_ens].isel(N_OBS=i_obs)

        for member in np.arange(N):
            #     print dsub['PSAL_LEVEL_QC'].sel(member=member).sel(N_OBS=i_obs).values
            #         plt.plot(x_mod.sel(member=member), y, 'r', linewidth=1)
            ax[iv].plot(x_mod.sel(N_MEMBER=member), y, '0.8', linewidth=1)
        hm, = ax[iv].plot(x_mod.mean('N_MEMBER'), y, 'k', label='Model Ensemble mean')
        ax[iv].plot(x_mod.mean('N_MEMBER') - x_mod.std('N_MEMBER'), y, 'k--')
        ax[iv].plot(x_mod.mean('N_MEMBER') + x_mod.std('N_MEMBER'), y, 'k--')
        ho, = ax[iv].plot(x_obs, y, 'r', label='Observation')
        # ax.legend([hm, ho], loc=4)
        ax[iv].set_xlabel(("%s [%s]")%(dsub[vname_obs].attrs['long_name'],dsub[vname_obs].attrs['units']))
        ax[iv].set_xlabel(("[%s]")%(dsub[vname_obs].attrs['units']))
        ax[iv].set_ylabel('Depth (m)')
        ax[iv].set_title(dsub[vname_obs].attrs['long_name'])

    ax[iv].invert_yaxis()
    plt.suptitle(("Measurement date: %s\nWMO: %s") % (ts.to_pydatetime().strftime("%Y/%m/%d %H:%M"), wmo))
    return fig

def ensemble_profiles_one(dsub, i_obs, vname = "POTM"):
    """Plot ensemble profiles superimposed with observation profile
        vname = "POTM"
        vname = "PSAL"
    """
    # WMO of the platform:
    wmo = np.asmatrix(dsub['STATION_IDENTIFIER'].isel(N_OBS=i_obs).values)[0, 0].strip()
    # Datetime of the measurement:
    jd = np.asmatrix(dsub['JULD'].isel(N_OBS=i_obs).values)[0, 0]
    ts = pd.Timestamp('1950-01-01') + pd.Timedelta(days=jd)
    # ts.to_pydatetime().strftime("%Y%m%d%H%M%S")

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
    plt.title(("Measurement date: %s\nWMO: %s") % (ts.to_pydatetime().strftime("%Y/%m/%d %H:%M"), wmo))
    return fig

def ensemble_profile_pdf(dsub, i_obs, N=5, vname="POTM", dpi=80, hspace=-.35, obscolor='axis'):
    """Plot the ensemble pdf of a vertical profile
        N is the number of depth levels
        vname = "POTM"
        vname = "PSAL"
        obscolor='axis'
        obscolor='r'
        hspace=-.35
        hspace=0
        dpi=80
        dpi=300
        N=10
        N=15
    """
    # WMO of the platform:
    wmo = np.asmatrix(dsub['STATION_IDENTIFIER'].isel(N_OBS=i_obs).values)[0, 0].strip()
    # Datetime of the measurement:
    jd = np.asmatrix(dsub['JULD'].isel(N_OBS=i_obs).values)[0, 0]
    ts = pd.Timestamp('1950-01-01') + pd.Timedelta(days=jd)
    # ts.to_pydatetime().strftime("%Y%m%d%H%M%S")

    print "Create Pandas dataframe:"
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

    print "Observation info:"
    print np.unique(d['STATION_IDENTIFIER'])
    print np.unique(d['STATION_TYPE'])
    print np.unique(d['LATITUDE'])
    print np.unique(d['LONGITUDE'])
    # print np.unique(d['JULD'])
    # print np.unique(d['JULD_REFERENCE'])[0]
    print pr.utils.read_ncjuld(d)

    # print "Now plot the PDF"
    vname_pdf = ("%s_Hx") % (vname)
    vname_obs = ("%s_OBS") % (vname)

    # plt.figure(dpi=300)
    sns.set(style="white", rc={"axes.facecolor": (0, 0, 0, 0)})
    # sns.set(style="dark", rc={"axes.facecolor": (0, 0, 0, 0)})

    # Initialize the FacetGrid object
    pal = sns.cubehelix_palette(N, rot=-.05, light=.7)
    g = sns.FacetGrid(df, row="DEPTH", hue="DEPTH", aspect=7, size=1, palette=pal, sharey=False)
    # g = sns.FacetGrid(df, row="DEPTH", hue="DEPTH", aspect=5, size=1, sharey=False)
    # g = sns.FacetGrid(df, row="DEPTH", hue="DEPTH", aspect=5, size=1, palette="Reds", sharey=False)
    # g = sns.FacetGrid(df, row="DEPTH", hue="DEPTH", aspect=5, size=1, palette="RdBu", sharey=False)

    # Set figure resolution:
    g.fig.set_dpi(dpi)

    # Draw the densities in a few steps
    g.map(sns.distplot, vname_pdf, hist=True, color="k")
    # g.map(sns.kdeplot, vname_pdf, clip_on=False, shade=True, alpha=1, lw=1.5, bw=.2)
    # g.map(sns.kdeplot, vname_pdf, clip_on=False, color="w", lw=2, bw=.2)
    # g.map(sns.kdeplot, vname_pdf, cumulative=True, clip_on=False, color="k", lw=2, bw=.2)
    # g.map(sns.distplot, vname_pdf, kde=True, color="b")

    # Define and use a simple function to label the plot in axes coordinates
    def label(x, color, label):
        ax = plt.gca()
        ax.text(0, .2, ("z=%0.1fm") % (float(label)), fontweight="bold", color=color,
                ha="left", va="center", transform=ax.transAxes)
    g.map(label, vname_pdf)

    # Add observations segments:
    def obsline(x, color, label):
        ax = plt.gca()
        if obscolor == 'axis':
            plt.axvline(np.unique(x), np.min(ax.get_ylim()), np.max(ax.get_ylim()), color=color)
        else:
            plt.axvline(np.unique(x), np.min(ax.get_ylim()), np.max(ax.get_ylim()), color=obscolor)
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
    g.fig.suptitle(("%s\nMeasurement: %s\nWMO: %s") % (dsub[vname_pdf].attrs['long_name'],
                                 ts.to_pydatetime().strftime("%Y/%m/%d %H:%M"), wmo))

    return g
