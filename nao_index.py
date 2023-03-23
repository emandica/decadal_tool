#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Feb 28 17:52:45 2023

@author: ema
"""
import pandas as pd

import os.path

import xarray as xr

from eofs.xarray import Eof

import numpy as np

import datetime

import matplotlib.pyplot as plt
import matplotlib.ticker as mticker
import matplotlib.path as mpath
plt.rcParams.update({'font.size': 22})

import cartopy.crs as ccrs

import data_agg as da
import level_selection as ls
import nao_bootstrap as nbs
import constants as c

def pcs_check(var_eofs,var_pcs):
    
    if var_eofs[0].sel(**{'lat':slice(63,70),'lon':slice(-25,-16)}).mean('lat').mean('lon') > 0:
        var_eofs=-var_eofs
        var_pcs=-var_pcs
        
    return var_eofs, var_pcs
        
    
def nao_index_eofs(var):
    var.coords['lon']=(var.coords['lon'] + 180) % 360 - 180
    var=var.sortby(var.lon)
    LatIndexer, LonIndexer = 'lat', 'lon'
    var = var.sel(**{LatIndexer: slice(20,80),
                     LonIndexer: slice(-90,40)})
    
    wgts = np.sqrt(np.cos(np.deg2rad(var.lat.values)))[...,np.newaxis]
    solver = Eof(var, weights=wgts)
    pcs = solver.pcs(npcs=1, pcscaling=1)
    eofs = solver.eofsAsCorrelation(neofs=1)
    variance_fractions = solver.varianceFraction(neigs=1)
    
    return pcs, eofs, variance_fractions

def nao_loop(dset):
    nao_pcs=[]
    nao_eofs=[]
    nao_var=[]
    for i in range(10):
        dset_mem = dset.isel({'member':i})
        anci1,anci2,anci3=nao_index_eofs(dset_mem)
        nao_pcs.append(anci1)
        nao_eofs.append(anci2)
        nao_var.append(anci3)
    nao_pcs = xr.concat(nao_pcs, dim='member')
    nao_eofs = xr.concat(nao_eofs, dim='member')
    nao_var = xr.concat(nao_var, dim = 'member')
    nao_pcs = nao_pcs.squeeze()
    nao_var = nao_var.squeeze()
    return nao_pcs, nao_eofs, nao_var

def index_plot(var,proj,levels,title):
    fig = plt.figure(figsize=[12,6])
    ax = fig.add_subplot(111,projection = proj)
    ax.coastlines()
    ax.set_global()
   
    p=var.plot.contourf(ax=ax,levels=levels,transform = ccrs.PlateCarree(), label=c.REF,add_labels=False,add_colorbar=False)
    
    # add separate colorbar
    cb = plt.colorbar(p, ticks=levels, shrink=1, extend='both')
    cb.ax.tick_params(labelsize=18)
    
    #Drow gridlines and adjust labels
    gl = p.axes.gridlines(crs=ccrs.PlateCarree(), draw_labels=False)
    gl.top_labels = False
    gl.right_labels = False

    ax.set_title(title, fontsize=16)

#%%

def NAO(lead_exp,season):
    if os.path.isfile(c.RUN_DIR+c.NAME_CTL+'_'+c.VAR+'_lead'+str(lead_exp)+'_'+season+'_s'+str(c.T_START)+'.nc') == False:
        da.aggr_datasets(lead_exp, season)
    if os.path.isfile(c.RUN_DIR+c.NAME_CTL+'_'+c.VAR+'_lead'+str(lead_exp)+'_'+season+'_s'+str(c.T_START)+'_level_'+str(c.PLEV)+'.nc') == False:            
        ls.level_sel(lead_exp, season)
        
    ctl = xr.open_dataset(c.RUN_DIR+c.NAME_CTL+'_'+c.VAR+'_lead'+str(lead_exp)+'_'+season+'_s'+str(c.T_START)+'_level_'+str(c.PLEV)+'.nc', chunks={'lon':'auto','lat':'auto'})
    sens = xr.open_dataset(c.RUN_DIR+c.NAME_SENS+'_'+c.VAR+'_lead'+str(lead_exp)+'_'+season+'_s'+str(c.T_START)+'_level_'+str(c.PLEV)+'.nc', chunks={'lon':'auto','lat':'auto'})
    era = xr.open_dataset(c.RUN_DIR+c.REF+'_'+c.VAR+'_lead'+str(lead_exp)+'_'+season+'_s'+str(c.T_START)+'_level_'+str(c.PLEV)+'.nc', chunks={'lon':'auto','lat':'auto'})
    era= era.sel(time=slice(ctl.time[0],ctl.time[-1]))
    
    #anomalies
    ctl = ctl - ctl.sel(time=slice(str(c.T_START+len(lead_exp)-1), ctl.time[-len(lead_exp)])).mean('time')
    sens = sens - sens.sel(time=slice(str(c.T_START+len(lead_exp)-1), ctl.time[-len(lead_exp)])).mean('time')
    era = era - era.sel(time=slice(str(c.T_START+len(lead_exp)-1), ctl.time[-len(lead_exp)])).mean('time')
    
    ctl_em = ctl.mean('member')
    sens_em = sens.mean('member')
#%% 
    nao_ctl_pcs, nao_ctl_eofs, nao_ctl_var=nao_loop(ctl[c.VAR])
    ctl_std_pcs = nao_ctl_pcs.std('member')
    ctl_std_eofs = nao_ctl_eofs.std('member')
    
    nao_sens_pcs, nao_sens_eofs,nao_sens_var=nao_loop(sens[c.VAR])
    sens_std_pcs = nao_sens_pcs.std('member')
    sens_std_eofs = nao_sens_eofs.std('member')

    nao_ctl_em_pcs, nao_ctl_em_eofs,nao_ctl_em_var = nao_index_eofs(ctl_em[c.VAR])
    nao_ctl_em_pcs = nao_ctl_em_pcs.squeeze()
    nao_ctl_em_eofs = nao_ctl_em_eofs
    
    nao_sens_em_pcs,nao_sens_em_eofs, nao_sens_em_var = nao_index_eofs(sens_em[c.VAR])
    nao_sens_em_pcs = nao_sens_em_pcs.squeeze()
    nao_sens_em_eofs = nao_sens_em_eofs
    
    nao_era_pcs,nao_era_eofs,nao_era_var = nao_index_eofs(era[c.VAR])
    
    #trass PCS sign
    nao_era_eofs,nao_era_pcs=pcs_check(nao_era_eofs,nao_era_pcs)
    nao_ctl_em_eofs,nao_ctl_em_pcs=pcs_check(nao_ctl_em_eofs,nao_ctl_em_pcs)
    nao_sens_em_eofs,nao_sens_em_pcs=pcs_check(nao_sens_em_eofs,nao_sens_em_pcs)
        
    
#%%save file
    nao_ctl_pcs.to_netcdf(c.RUN_DIR+c.NAME_CTL+'_nao_index_lead'+str(lead_exp)+'_'+season+'.nc')
    nao_ctl_eofs.to_netcdf(c.RUN_DIR+c.NAME_CTL+'_nao_eof1_lead'+str(lead_exp)+'_'+season+'.nc')
    nao_ctl_var.to_netcdf(c.RUN_DIR+c.NAME_CTL+'_variance_lead'+str(lead_exp)+'_'+season+'.nc')

    nao_ctl_em_pcs.to_netcdf(c.RUN_DIR+c.NAME_CTL+'_em_nao_index_lead'+str(lead_exp)+'_'+season+'.nc')
    nao_ctl_em_eofs.to_netcdf(c.RUN_DIR+c.NAME_CTL+'_em_nao_eof1_lead'+str(lead_exp)+'_'+season+'.nc')
    nao_ctl_em_var.to_netcdf(c.RUN_DIR+c.NAME_CTL+'_em_variance_lead'+str(lead_exp)+'_'+season+'.nc')

    nao_sens_pcs.to_netcdf(c.RUN_DIR+c.NAME_SENS+'_nao_index_lead'+str(lead_exp)+'_'+season+'.nc')
    nao_sens_eofs.to_netcdf(c.RUN_DIR+c.NAME_SENS+'_nao_eof1_lead'+str(lead_exp)+'_'+season+'.nc')
    nao_sens_var.to_netcdf(c.RUN_DIR+c.NAME_SENS+'_variance_lead'+str(lead_exp)+'_'+season+'.nc')

    nao_sens_em_pcs.to_netcdf(c.RUN_DIR+c.NAME_SENS+'_em_nao_index_lead'+str(lead_exp)+'_'+season+'.nc')
    nao_sens_em_eofs.to_netcdf(c.RUN_DIR+c.NAME_SENS+'_em_nao_eof1_lead'+str(lead_exp)+'_'+season+'.nc')
    nao_sens_em_var.to_netcdf(c.RUN_DIR+c.NAME_SENS+'_em_variance_lead'+str(lead_exp)+'_'+season+'.nc')
    
    nao_era_pcs.to_netcdf(c.RUN_DIR+c.REF+'_nao_index_lead'+str(lead_exp)+'_'+season+'.nc')
    nao_era_eofs.to_netcdf(c.RUN_DIR+c.REF+'_nao_eof1_lead'+str(lead_exp)+'_'+season+'.nc')
    nao_era_var.to_netcdf(c.RUN_DIR+c.REF+'_variance_lead'+str(lead_exp)+'_'+season+'.nc')

    ###
    ctl_nao_cor = xr.corr(nao_ctl_em_pcs, nao_era_pcs, dim='time').to_numpy().round(2)
    sens_nao_cor = xr.corr(nao_sens_em_pcs, nao_era_pcs, dim='time').to_numpy().round(2)
    
    d_corr = sens_nao_cor - ctl_nao_cor
    
    quant_corr = nbs.NAO_bootstrap(lead_exp, season)
    
    if d_corr < quant_corr[0] or d_corr>quant_corr[-1]:
        sign = True
    else:
        sign = False
#%%
    fig = plt.figure(figsize=[12,8])
    ax = fig.add_subplot(111)
    
    nao_ctl_em_pcs.plot.line(x='time',c='red',label=c.NAME_CTL+' r='+str(ctl_nao_cor))
    _ = ax.fill_between(nao_ctl_em_pcs.time,nao_ctl_em_pcs-ctl_std_pcs,nao_ctl_em_pcs+ctl_std_pcs,alpha=0.2,color='red') 

    nao_sens_em_pcs.plot(x='time',c='blue',label=c.NAME_SENS+' r='+str(sens_nao_cor))
    _ = ax.fill_between(nao_sens_em_pcs.time,nao_sens_em_pcs-sens_std_pcs,nao_sens_em_pcs+sens_std_pcs,alpha=0.2,color='blue') 
    
    nao_era_pcs.plot(label=c.REF,c='k')

    ax.set_ylabel=('NAO index')
    ax.set_title('NAO index, lead '+str(lead_exp)+', '+season+', corr_sign=95% '+str(sign) )
    
    ax.set_ylim([-4,4])
    
    plt.grid()
    
    ax.legend(frameon=False)
    plt.savefig(c.OUT_DIR + 'NAO_index_lead_'+str(lead_exp) + '_' + season + c.REF + '.jpg', dpi=300, bbox_inches='tight')    

#%%
    levels= [-1,-0.8,-0.6,-0.4,-0.2,0.2,0.4,0.6,0.8,1]
    
#%%    
    proj = ccrs.Orthographic(central_longitude=-20, central_latitude=60)

    title = 'NAO_EOF1_lead_'+str(lead_exp) + '_' + season +'_'+ c.NAME_CTL+'_var='+str((nao_ctl_var[0].values*100).round(1))
    index_plot(nao_ctl_em_eofs[0], proj, levels, title)
    plt.savefig(c.OUT_DIR + '_'+c.NAME_CTL+'_NAO_EOF1_lead_'+str(lead_exp) + '_' + season +'_'+ c.REF + '.jpg', dpi=300, bbox_inches='tight')    
    
#%%
    title='NAO_EOF1_lead_'+str(lead_exp) + '_' + season +'_'+ c.NAME_SENS+'_var='+str((nao_sens_var[0].values*100).round(1))
    index_plot(nao_sens_em_eofs[0], proj, levels, title)
    plt.savefig(c.OUT_DIR + '_'+c.NAME_SENS+'_NAO_EOF1_lead_'+str(lead_exp) + '_' + season +'_'+ c.REF + '.jpg', dpi=300, bbox_inches='tight')    
#%%
    title= 'NAO_EOF1_lead_'+str(lead_exp) + '_' + season +'_'+ c.REF+'_var='+str((nao_era_var[0].values*100).round(1))
    index_plot(nao_era_eofs[0], proj, levels, title)
    plt.savefig(c.OUT_DIR + '_'+c.REF+'_NAO_EOF1_lead_'+str(lead_exp) + '_' + season +'_'+ c.REF + '.jpg', dpi=300, bbox_inches='tight')    