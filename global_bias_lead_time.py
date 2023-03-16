#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Feb 27 18:05:38 2023

@author: ema
"""

import xarray as xr
import os.path

import numpy as np
import matplotlib.pyplot as plt
plt.rcParams.update({'font.size': 22})
import constants as c
import data_agg as da
import level_selection as ls

def global_bias(season):
    for i in [[0],[1],[2],[3],[4]]:
        if os.path.isfile(c.RUN_DIR+c.NAME_CTL+'_'+c.VAR+'_lead'+str(i)+'_'+season+'_s'+str(c.T_START)+'.nc') == False:
            da.aggr_datasets(i, season)
        if os.path.isfile(c.RUN_DIR+c.NAME_CTL+'_'+c.VAR+'_lead'+str(i)+'_'+season+'_s'+str(c.T_START)+'_level_'+str(c.PLEV)+'.nc') == False:            
            ls.level_sel(i, season)
#%%    
    ctl_path_0 = c.RUN_DIR+c.NAME_CTL+'_'+c.VAR+'_lead[0]_'+season+'_s'+str(c.T_START)+'_level_'+str(c.PLEV)+'.nc'
    ctl_0 = xr.open_dataset(ctl_path_0, chunks={'lon':'auto','lat':'auto'})
    ctl_0=ctl_0.mean('time')
    
    ctl_path_1 = c.RUN_DIR+c.NAME_CTL+'_'+c.VAR+'_lead[1]_'+season+'_s'+str(c.T_START)+'_level_'+str(c.PLEV)+'.nc'
    ctl_1 = xr.open_dataset(ctl_path_1, chunks={'lon':'auto','lat':'auto'})
    ctl_1=ctl_1.mean('time')
    
    ctl_path_2 = c.RUN_DIR+c.NAME_CTL+'_'+c.VAR+'_lead[2]_'+season+'_s'+str(c.T_START)+'_level_'+str(c.PLEV)+'.nc'
    ctl_2 = xr.open_dataset(ctl_path_2, chunks={'lon':'auto','lat':'auto'})
    ctl_2=ctl_2.mean('time')
    
    ctl_path_3 = c.RUN_DIR+c.NAME_CTL+'_'+c.VAR+'_lead[3]_'+season+'_s'+str(c.T_START)+'_level_'+str(c.PLEV)+'.nc'
    ctl_3 = xr.open_dataset(ctl_path_3, chunks={'lon':'auto','lat':'auto'})
    ctl_3=ctl_3.mean('time')
    
    ctl_path_4 = c.RUN_DIR+c.NAME_CTL+'_'+c.VAR+'_lead[4]_'+season+'_s'+str(c.T_START)+'_level_'+str(c.PLEV)+'.nc'
    ctl_4 = xr.open_dataset(ctl_path_4, chunks={'lon':'auto','lat':'auto'})
    ctl_4=ctl_4.mean('time')
#%%
    sens_path_0 = c.RUN_DIR+c.NAME_SENS+'_'+c.VAR+'_lead[0]_'+season+'_s'+str(c.T_START)+'_level_'+str(c.PLEV)+'.nc'
    sens_0 = xr.open_dataset(sens_path_0, chunks={'lon':'auto','lat':'auto'})
    sens_0=sens_0.mean('time')
    
    sens_path_1 = c.RUN_DIR+c.NAME_SENS+'_'+c.VAR+'_lead[1]_'+season+'_s'+str(c.T_START)+'_level_'+str(c.PLEV)+'.nc'
    sens_1 = xr.open_dataset(sens_path_1, chunks={'lon':'auto','lat':'auto'})
    sens_1=sens_1.mean('time')
    
    sens_path_2 = c.RUN_DIR+c.NAME_SENS+'_'+c.VAR+'_lead[2]_'+season+'_s'+str(c.T_START)+'_level_'+str(c.PLEV)+'.nc'
    sens_2 = xr.open_dataset(sens_path_2, chunks={'lon':'auto','lat':'auto'})
    sens_2=sens_2.mean('time')
    
    sens_path_3 = c.RUN_DIR+c.NAME_SENS+'_'+c.VAR+'_lead[3]_'+season+'_s'+str(c.T_START)+'_level_'+str(c.PLEV)+'.nc'
    sens_3 = xr.open_dataset(sens_path_3, chunks={'lon':'auto','lat':'auto'})
    sens_3=sens_3.mean('time')
    
    sens_path_4 = c.RUN_DIR+c.NAME_SENS+'_'+c.VAR+'_lead[4]_'+season+'_s'+str(c.T_START)+'_level_'+str(c.PLEV)+'.nc'
    sens_4 = xr.open_dataset(sens_path_4, chunks={'lon':'auto','lat':'auto'})    
    sens_4=sens_4.mean('time')

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
#%%
    a52o_land_file = '/mnt/d/00-dataset/02-confess/03-confess_data/a52o/land/fc00_r.nc'
    dset = xr.open_dataset(a52o_land_file)
    dset=dset['LSM']
    
    ctl = ctl[c.VAR]
    sens = sens[c.VAR]
    era = era[c.VAR]
    
    for mask in ['global','land','ocean']:
    
        if mask=='land':
            nctl = xr.where(dset.mean('time')>0.5, ctl, np.nan)
            nsens = xr.where(dset.mean('time')>0.5, sens, np.nan)
            nera = xr.where(dset.mean('time')>0.5, era, np.nan)
        elif mask=='ocean':
            nctl = xr.where(dset.mean('time')<0.5, ctl, np.nan)
            nsens = xr.where(dset.mean('time')<0.5, sens, np.nan)
            nera = xr.where(dset.mean('time')<0.5, era, np.nan)
        else:
            nctl=ctl
            nsens=sens
            nera=era

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

        nctl = nctl.sel(**{LatIndexer: slice(lat_min,lat_max) , LonIndexer: slice(lon_min,lon_max)}).weighted(weights).mean(LatIndexer).mean(LonIndexer)
        
        nsens = nsens.sel(**{LatIndexer: slice(lat_min,lat_max),LonIndexer:slice(lon_min,lon_max)}).weighted(weights).mean(LatIndexer).mean(LonIndexer)
    
        nera = nera.sel(**{LatIndexer: slice(lat_min,lat_max),LonIndexer:slice(lon_min,lon_max)}).weighted(weights).mean(LatIndexer).mean(LonIndexer)
    
        nctl = nctl-nera
        nsens = nsens-nera   

        nctl = nctl.mean('time')
        nsens = nsens.mean('time')
    
        ctl_std=nctl.std('member')
        nctl=nctl.mean('member')

        sens_std=nsens.std('member')    
        nsens=nsens.mean('member')
    
        fig = plt.figure(figsize=[12,8])
        ax = fig.add_subplot(111)
        
        p = nctl.plot.line(x='lead_year',label=c.NAME_CTL)
        _ = nsens.plot(label=c.NAME_SENS)
        _ = ax.fill_between(nctl.lead_year,nctl-ctl_std,nctl+ctl_std,alpha=0.5) 
        _ = ax.fill_between(nsens.lead_year,nsens-sens_std,nsens+sens_std,alpha=0.5) 
        plt.legend()
    
        ax.set_xticks([0,1,2,3,4])
        ax.set_xlim([0,4])
        ax.set_ylabel(c.UNITS)
        #if c.VAR == 'ua':
            #ax.set_ylim([-0.2,0])
        ax.set_title(c.NAME_SENS+'-'+c.NAME_CTL+'_'+mask+'_bias_'+c.VAR+'('+c.UNITS+')_' + season + '_s'+str(c.T_START)+ '_' + c.REF)
    
        plt.savefig(c.OUT_DIR + c.NAME_SENS+'-'+c.NAME_CTL+'_'+mask+'_BIAS_'+c.VAR+ '_' + season + '_s'+str(c.T_START)+'_' + c.REF +'.jpg', dpi=300, bbox_inches='tight')