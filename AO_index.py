#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Mar 14 11:07:04 2023

@author: ema
"""
import pandas as pd

import xarray as xr

from eofs.xarray import Eof

import matplotlib.pyplot as plt
plt.rcParams.update({'font.size': 22})

import numpy as np

import datetime

import constants as c

def ao_index_eofs(var):
    var.coords['lon']=(var.coords['lon'] + 180) % 360 - 180
    var=var.sortby(var.lon)
    LatIndexer, LonIndexer = 'lat', 'lon'
    var = var.sel(**{LatIndexer: slice(20,90)})
    var = var- var.mean('time')
    
    coslat = np.cos(var.coords['lat'].values).clip(0., 1.)
    wgts = np.sqrt(coslat)[..., np.newaxis]
    solver = Eof(var, weights=wgts)
    eofs = solver.pcs(npcs=1, pcscaling=1)
    #eofs = solver.eofs(neofs=1)
    
    return eofs

def ao_loop(dset):
    ao=[]
    for i in range(10):
        dset_mem = dset.isel({'member':i})
        anci=ao_index_eofs(dset_mem)
        ao.append(anci)
    ao = xr.concat(ao, dim='member')
    ao = ao.squeeze()
    return ao

def AO(lead_exp,season):
    ctl = xr.open_dataset(c.RUN_DIR+c.NAME_CTL+'_'+c.VAR+'_lead'+str(lead_exp)+'_'+season+'_s'+str(c.T_START)+'_level_'+str(c.PLEV)+'.nc', chunks={'lon':'auto','lat':'auto'})
    sens = xr.open_dataset(c.RUN_DIR+c.NAME_SENS+'_'+c.VAR+'_lead'+str(lead_exp)+'_'+season+'_s'+str(c.T_START)+'_level_'+str(c.PLEV)+'.nc', chunks={'lon':'auto','lat':'auto'})
    era = xr.open_dataset(c.RUN_DIR+c.REF+'_'+c.VAR+'_lead'+str(lead_exp)+'_'+season+'_s'+str(c.T_START)+'_level_'+str(c.PLEV)+'.nc', chunks={'lon':'auto','lat':'auto'})
    era= era.sel(time=slice(ctl.time[0],ctl.time[-1]))
    
    ctl_em = ctl.mean('member')
    sens_em = sens.mean('member')
    
#%%
    ao_ctl=ao_loop(ctl[c.VAR])
    ctl_std = ao_ctl.std('member')
    
    ao_sens=ao_loop(sens[c.VAR])
    sens_std = ao_sens.std('member')
    
    ao_ctl_em = ao_index_eofs(ctl_em[c.VAR]).squeeze() 
    ao_sens_em = ao_index_eofs(sens_em[c.VAR]).squeeze() 
    ao_era = ao_index_eofs(era[c.VAR])
    
    ctl_ao_cor = xr.corr(ao_ctl_em, ao_era, dim='time').to_numpy().round(2)
    sens_ao_cor = xr.corr(ao_sens_em, ao_era, dim='time').to_numpy().round(2)
    
#%%
    fig = plt.figure(figsize=[12,6])
    ax = fig.add_subplot(111)

    ao_ctl_em.plot.line(x='time',c='red',label=c.NAME_CTL+' r='+str(ctl_ao_cor))
    _ = ax.fill_between(ao_ctl_em.time,ao_ctl_em-ctl_std,ao_ctl_em+ctl_std,alpha=0.2,color='red') 

    ao_sens_em.plot(x='time',c='blue',label=c.NAME_SENS+' r='+str(sens_ao_cor))
    _ = ax.fill_between(ao_sens_em.time,ao_sens_em-sens_std,ao_sens_em+sens_std,alpha=0.2,color='blue') 
        
    ao_era.plot(label=c.REF,c='k')

    ax.set_ylabel=('AO index')
    ax.set_title('AO index, lead '+str(lead_exp)+', '+season)
        
    ax.set_ylim([-2.5,2.5])
    
    plt.legend()
    plt.savefig(c.OUT_DIR + 'AO_index_lead_'+str(lead_exp) + '_' + season + c.REF + '.jpg', dpi=300, bbox_inches='tight')    
    