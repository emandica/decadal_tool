#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jan 26 11:17:14 2023

@author: ema
"""
import xarray as xr

import constants as c
import open_dataset as od
import generic_routines as cr

###############################################################################
def aggr_datasets(lead_exp, season):
#%%
    """
    open datasets
    """
    
    ctl, sens, ref = od.open_dataset_and_seasonal_selection(lead_exp, season)
    
    """
    bias correction
    """
    ctl = ctl-ctl.mean('member').mean('time')
    sens = sens-sens.mean('member').mean('time')

    """
    lead mean
    """
    ctl = ctl.mean('lead')
    sens = sens.mean('lead')
    ref = ref.rolling(time=len(lead_exp)).mean()
    ref = ref.shift(time=-round(len(lead_exp)-1))
    
#%%
    """
    save dataset
    """
    ctl.to_netcdf(c.RUN_DIR+c.NAME_CTL+'_'+c.VAR+'_lead'+str(lead_exp)+'_'+season+'_s'+str(c.T_START)+'.nc')
    sens.to_netcdf(c.RUN_DIR+c.NAME_SENS+'_'+c.VAR+'_lead'+str(lead_exp)+'_'+season+'_s'+str(c.T_START)+'.nc')
    ref.to_netcdf(c.RUN_DIR+c.REF+'_'+c.VAR+'_lead'+str(lead_exp)+'_'+season+'_s'+str(c.T_START)+'.nc')

#%%
    """
    anomaly calculation
    """
    #ctl = ctl - ctl.sel(time=slice(str(c.T_START+len(lead_exp)-1), ctl.time[-len(lead_exp)])).mean('member').mean('time')
    ctl = ctl - ctl.sel(time=slice(str(c.T_START+len(lead_exp)-1), ctl.time[-len(lead_exp)])).mean('time')
    #sens = sens - sens.sel(time=slice(str(c.T_START+len(lead_exp)-1), ctl.time[-len(lead_exp)])).mean('member').mean('time')
    sens = sens - sens.sel(time=slice(str(c.T_START+len(lead_exp)-1), ctl.time[-len(lead_exp)])).mean('time')
    ref = ref - ref.sel(time=slice(str(c.T_START+len(lead_exp)-1), ctl.time[-len(lead_exp)])).mean('time')

#%%
    """
    save dataset
    """
    ctl.to_netcdf(c.RUN_DIR+c.NAME_CTL+'_'+c.VAR+'_lead'+str(lead_exp)+'_'+season+'_s'+str(c.T_START)+'_anomaly.nc')
    sens.to_netcdf(c.RUN_DIR+c.NAME_SENS+'_'+c.VAR+'_lead'+str(lead_exp)+'_'+season+'_s'+str(c.T_START)+'_anomaly.nc')
    ref.to_netcdf(c.RUN_DIR+c.REF+'_'+c.VAR+'_lead'+str(lead_exp)+'_'+season+'_s'+str(c.T_START)+'_anomaly.nc')
    