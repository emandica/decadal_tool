#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Feb 27 21:47:44 2023

@author: ema
"""

import xarray as xr

import numpy as np
import matplotlib.pyplot as plt
plt.rcParams.update({'font.size': 22})

import constants as c

def global_correlation(season):
    
#%%    
    ctl_path_0 = c.RUN_DIR+c.NAME_CTL+'_'+c.VAR+'_lead[0]_'+season+'_s'+str(c.T_START)+'_level_'+str(c.PLEV)+'.nc'
    ctl_0 = xr.open_dataset(ctl_path_0, chunks={'lon':'auto','lat':'auto'})
    
    ctl_path_1 = c.RUN_DIR+c.NAME_CTL+'_'+c.VAR+'_lead[1]_'+season+'_s'+str(c.T_START)+'_level_'+str(c.PLEV)+'.nc'
    ctl_1 = xr.open_dataset(ctl_path_1, chunks={'lon':'auto','lat':'auto'})
    
    ctl_path_2 = c.RUN_DIR+c.NAME_CTL+'_'+c.VAR+'_lead[2]_'+season+'_s'+str(c.T_START)+'_level_'+str(c.PLEV)+'.nc'
    ctl_2 = xr.open_dataset(ctl_path_2, chunks={'lon':'auto','lat':'auto'})
    
    ctl_path_3 = c.RUN_DIR+c.NAME_CTL+'_'+c.VAR+'_lead[3]_'+season+'_s'+str(c.T_START)+'_level_'+str(c.PLEV)+'.nc'
    ctl_3 = xr.open_dataset(ctl_path_3, chunks={'lon':'auto','lat':'auto'})
    
    ctl_path_4 = c.RUN_DIR+c.NAME_CTL+'_'+c.VAR+'_lead[4]_'+season+'_s'+str(c.T_START)+'_level_'+str(c.PLEV)+'.nc'
    ctl_4 = xr.open_dataset(ctl_path_4, chunks={'lon':'auto','lat':'auto'})
#%%
    sens_path_0 = c.RUN_DIR+c.NAME_SENS+'_'+c.VAR+'_lead[0]_'+season+'_s'+str(c.T_START)+'_level_'+str(c.PLEV)+'.nc'
    sens_0 = xr.open_dataset(sens_path_0, chunks={'lon':'auto','lat':'auto'})
    
    sens_path_1 = c.RUN_DIR+c.NAME_SENS+'_'+c.VAR+'_lead[1]_'+season+'_s'+str(c.T_START)+'_level_'+str(c.PLEV)+'.nc'
    sens_1 = xr.open_dataset(sens_path_1, chunks={'lon':'auto','lat':'auto'})
    
    sens_path_2 = c.RUN_DIR+c.NAME_SENS+'_'+c.VAR+'_lead[2]_'+season+'_s'+str(c.T_START)+'_level_'+str(c.PLEV)+'.nc'
    sens_2 = xr.open_dataset(sens_path_2, chunks={'lon':'auto','lat':'auto'})
    
    sens_path_3 = c.RUN_DIR+c.NAME_SENS+'_'+c.VAR+'_lead[3]_'+season+'_s'+str(c.T_START)+'_level_'+str(c.PLEV)+'.nc'
    sens_3 = xr.open_dataset(sens_path_3, chunks={'lon':'auto','lat':'auto'})
    
    sens_path_4 = c.RUN_DIR+c.NAME_SENS+'_'+c.VAR+'_lead[4]_'+season+'_s'+str(c.T_START)+'_level_'+str(c.PLEV)+'.nc'
    sens_4 = xr.open_dataset(sens_path_4, chunks={'lon':'auto','lat':'auto'})    

#%%
    era_path_0 =c.RUN_DIR+c.REF+'_'+c.VAR+'_lead[0]'+'_'+season+'_s'+str(c.T_START)+'_level_'+str(c.PLEV)+'.nc'
    era_0 = xr.open_dataset(era_path_0,chunks={'lon':'auto','lat':'auto'})
    era_path_1 =c.RUN_DIR+c.REF+'_'+c.VAR+'_lead[1]'+'_'+season+'_s'+str(c.T_START)+'_level_'+str(c.PLEV)+'.nc'
    era_1 = xr.open_dataset(era_path_1,chunks={'lon':'auto','lat':'auto'})
    era_path_2 =c.RUN_DIR+c.REF+'_'+c.VAR+'_lead[2]'+'_'+season+'_s'+str(c.T_START)+'_level_'+str(c.PLEV)+'.nc'
    era_2 = xr.open_dataset(era_path_2,chunks={'lon':'auto','lat':'auto'})
    era_path_3 =c.RUN_DIR+c.REF+'_'+c.VAR+'_lead[3]'+'_'+season+'_s'+str(c.T_START)+'_level_'+str(c.PLEV)+'.nc'
    era_3 = xr.open_dataset(era_path_3,chunks={'lon':'auto','lat':'auto'})
    era_path_4 =c.RUN_DIR+c.REF+'_'+c.VAR+'_lead[4]'+'_'+season+'_s'+str(c.T_START)+'_level_'+str(c.PLEV)+'.nc'
    era_4 = xr.open_dataset(era_path_4,chunks={'lon':'auto','lat':'auto'})

    ctl = xr.concat([ctl_0,ctl_1,ctl_2,ctl_3,ctl_4],dim = 'lead_year') 
    sens = xr.concat([sens_0,sens_1,sens_2,sens_3,sens_4],dim = 'lead_year') 
    era = xr.concat([era_0,era_1,era_2,era_3,era_4],dim = 'lead_year') 

    ctl=ctl-ctl.mean('time')
    sens=sens-sens.mean('time')
    era= era-era.mean('time')
#%%
        a52o_land_file = '/mnt/d/00-dataset/02-confess/03-confess_data/a52o/land/fc00_r.nc'
        dset = xr.open_dataset(a52o_land_file)
        dset=dset['LSM']
        
        if mask=='land':
            ctl = xr.where(dset.mean('time')>0.5, ctl, np.nan)
            sens = xr.where(dset.mean('time')>0.5, sens, np.nan)
            era = xr.where(dset.mean('time')>0.5, era, np.nan)
        elif mask=='ocean':
            ctl = xr.where(dset.mean('time')<0.5, ctl, np.nan)
            sens = xr.where(dset.mean('time')<0.5, sens, np.nan)
            era = xr.where(dset.mean('time')<0.5, era, np.nan)
            
#%%
    m_ctl_cor = xr.corr(ctl[c.VAR].mean('member'),era[c.VAR],dim='time')
    ctl_cor = xr.corr(ctl[c.VAR],era[c.VAR],dim='time')            
    
    m_sens_cor = xr.corr(sens[c.VAR].mean('member'),era[c.VAR],dim='time')
    sens_cor = xr.corr(sens[c.VAR],era[c.VAR],dim='time')            

#%%
    """
    global mean
    """
    lat_min = -90
    lat_max = 90

    lon_min = 0
    lon_max = 360


    LatIndexer, LonIndexer = 'lat', 'lon'
    weights = np.cos(np.deg2rad(ctl_0[c.VAR].lat))
    weights.name = "weights"

    ctl_cor = ctl_cor.sel(**{LatIndexer: slice(lat_min,lat_max) , LonIndexer: slice(lon_min,lon_max)}).weighted(weights).mean(LatIndexer).mean(LonIndexer)
    m_ctl_cor = m_ctl_cor.sel(**{LatIndexer: slice(lat_min,lat_max) , LonIndexer: slice(lon_min,lon_max)}).weighted(weights).mean(LatIndexer).mean(LonIndexer)
        
    sens_cor = sens_cor.sel(**{LatIndexer: slice(lat_min,lat_max),LonIndexer:slice(lon_min,lon_max)}).weighted(weights).mean(LatIndexer).mean(LonIndexer)
    m_sens_cor = m_sens_cor.sel(**{LatIndexer: slice(lat_min,lat_max),LonIndexer:slice(lon_min,lon_max)}).weighted(weights).mean(LatIndexer).mean(LonIndexer)
        
#%%
    ctl_std=ctl_cor.std('member')
    sens_std=sens_cor.std('member')
    
#%%
    fig = plt.figure(figsize=[12,8])
    ax = fig.add_subplot(111)
        
    p = m_ctl_cor.plot(ax=ax,label=c.NAME_CTL)
    _ = m_sens_cor.plot(label=c.NAME_SENS)
    _ = ax.fill_between(m_ctl_cor.lead_year,m_ctl_cor-ctl_std,m_ctl_cor+ctl_std,alpha=0.5) 
    _ = ax.fill_between(m_sens_cor.lead_year,m_sens_cor-sens_std,m_sens_cor+sens_std,alpha=0.5) 
    plt.legend()    
    
    ax.set_xticks([0,1,2,3,4])
    ax.set_xlim([0,4])
    ax.set_ylim([0,0.5])
    ax.set_title(c.NAME_SENS+'-'+c.NAME_CTL+'_'+mask+'_ACC_'+c.VAR+'('+c.UNITS+')_' + season + '_s'+str(c.T_START)+ '_' + c.REF)
    
    plt.savefig(c.OUT_DIR + c.NAME_SENS+'-'+c.NAME_CTL+'_'+mask+'_ACC_'+c.VAR+ '_' + season + '_s'+str(c.T_START)+'_' + c.REF +'.jpg', dpi=300, bbox_inches='tight')