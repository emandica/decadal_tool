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
plt.rcParams.update({'font.size': 22})

import cartopy.crs as ccrs

import data_agg as da
import level_selection as ls

import constants as c

def nao_index_eofs(var):
    var.coords['lon']=(var.coords['lon'] + 180) % 360 - 180
    var=var.sortby(var.lon)
    LatIndexer, LonIndexer = 'lat', 'lon'
    var = var.sel(**{LatIndexer: slice(20,80),
                     LonIndexer: slice(-90,40)})
    var = var- var.mean('time')
    
    wgts = np.sqrt(np.cos(np.deg2rad(var.lat.values)))[...,np.newaxis]
    solver = Eof(var, weights=wgts)
    pcs = solver.pcs(npcs=1, pcscaling=1)
    eofs = solver.eofsAsCorrelation(neofs=1)
    
    return pcs, eofs

def nao_index_box(var):
    var.coords['lon']=(var.coords['lon'] + 180) % 360 - 180
    var=var.sortby(var.lon)    
    LatIndexer, LonIndexer = 'lat', 'lon'
    
    var = var - var.mean('time')
    
    weights = np.cos(np.deg2rad(var.lat))
    weights.name = "weights"
    
    nao = var.sel(**{LatIndexer: slice(36,40),LonIndexer: slice(-28,-20)}).weighted(weights).mean(['lat','lon'])-var.sel(**{LatIndexer: slice(63,70),LonIndexer: slice(-25,-16)}).weighted(weights).mean(['lat','lon'])
    
    return nao

def nao_loop(dset):
    nao_pcs=[]
    nao_eofs=[]
    for i in range(10):
        dset_mem = dset.isel({'member':i})
        anci1,anci2=nao_index_eofs(dset_mem)
        nao_pcs.append(anci1)
        nao_eofs.append(anci2)
    nao_pcs = xr.concat(nao_pcs, dim='member')
    nao_eofs = xr.concat(nao_eofs, dim='member')
    nao_pcs = nao_pcs.squeeze()
    nao_eofs = nao_eofs
    return nao_pcs, nao_eofs

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
    
    ctl_em = ctl.mean('member')
    sens_em = sens.mean('member')
#%% 
    nao_ctl_pcs, nao_ctl_eofs=nao_loop(ctl[c.VAR])
    ctl_std_pcs = nao_ctl_pcs.std('member')
    ctl_std_eofs = nao_ctl_eofs.std('member')
    
    nao_sens_pcs, nao_sens_eofs=nao_loop(sens[c.VAR])
    sens_std_pcs = nao_sens_pcs.std('member')
    sens_std_eofs = nao_sens_eofs.std('member')

    nao_ctl_em_pcs, nao_ctl_em_eofs = nao_index_eofs(ctl_em[c.VAR])
    nao_ctl_em_pcs = nao_ctl_em_pcs.squeeze()
    nao_ctl_em_eofs = nao_ctl_em_eofs
    
    nao_sens_em_pcs,nao_sens_em_eofs = nao_index_eofs(sens_em[c.VAR])
    nao_sens_em_pcs = nao_sens_em_pcs.squeeze()
    nao_sens_em_eofs = nao_sens_em_eofs
    
    nao_era_pcs,nao_era_eofs = nao_index_eofs(era[c.VAR])

    ctl_nao_cor = xr.corr(nao_ctl_em_pcs, nao_era_pcs, dim='time').to_numpy().round(2)
    sens_nao_cor = xr.corr(nao_sens_em_pcs, nao_era_pcs, dim='time').to_numpy().round(2)
#%%
    fig = plt.figure(figsize=[12,6])
    ax = fig.add_subplot(111)
    
    nao_ctl_em_pcs.plot.line(x='time',c='red',label=c.NAME_CTL+' r='+str(ctl_nao_cor))
    _ = ax.fill_between(nao_ctl_em_pcs.time,nao_ctl_em_pcs-ctl_std_pcs,nao_ctl_em_pcs+ctl_std_pcs,alpha=0.2,color='red') 

    nao_sens_em_pcs.plot(x='time',c='blue',label=c.NAME_SENS+' r='+str(sens_nao_cor))
    _ = ax.fill_between(nao_sens_em_pcs.time,nao_sens_em_pcs-sens_std_pcs,nao_sens_em_pcs+sens_std_pcs,alpha=0.2,color='blue') 
    
    nao_era_pcs.plot(label=c.REF,c='k')

    ax.set_ylabel=('NAO index')
    ax.set_title('NAO index, lead '+str(lead_exp)+', '+season)
    
    ax.set_ylim([-2.5,2.5])
    
    plt.legend()
    plt.savefig(c.OUT_DIR + 'NAO_index_lead_'+str(lead_exp) + '_' + season + c.REF + '.jpg', dpi=300, bbox_inches='tight')    

#%%
    proj = ccrs.Orthographic(central_longitude=-20, central_latitude=60)
    fig = plt.figure(figsize=[12,6])
    ax = fig.add_subplot(111,projection = proj)
    ax.coastlines()
    ax.set_global()
    
    nao_ctl_em_eofs[0].plot.contourf(ax=ax,transform = ccrs.PlateCarree(), label=c.NAME_CTL+' r='+str(ctl_nao_cor))

    #nao_sens_em_eofs.plot(label=c.NAME_SENS+' r='+str(sens_nao_cor))
    
    #nao_era_pcs.plot(label=c.REF,c='k')

    #ax.set_ylabel=('NAO index')
    #ax.set_title('NAO index, lead '+str(lead_exp)+', '+season)
    
    #ax.set_ylim([-2.5,2.5])
    ax.set_title('EOF1 expressed as covariance, lead '+str(lead_exp), fontsize=16)
    #plt.legend()
    plt.savefig(c.OUT_DIR + 'NAO_EOF1_lead_'+str(lead_exp) + '_' + season + c.REF + '.jpg', dpi=300, bbox_inches='tight')    
