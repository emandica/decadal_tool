#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jan 18 11:47:05 2023

@author: ema
"""
import xarray as xr
import generic_routines as ncr
import constants as c

def open_dataset_and_seasonal_selection(lead_exp, season):
    
    """
Load datasets
    """    
    ctl = ncr.lead_aggregation(lead_exp, c.VAR, c.CTL_NUMBER, c.NAME_CTL,
                               c.T_START,c.T_END,
                               no_member=c.no_member,
                               )
    print('aggregation ctl: OK')
    
    sens = ncr.lead_aggregation(lead_exp, c.VAR, c.SENS_NUMBER, c.NAME_SENS,
                                c.T_START,c.T_END,
                                no_member=c.no_member,
                                )
    print('aggregation sens: OK')

    if c.REF=='era':
        ref = ncr.open_era5(c.DIRECTORY_ERA, c.FILE_ERA, c.VAR)
        ref = ref.sel(time=slice(ctl.time[0].dt.strftime("%Y-%m"),
                                           ref.time[-1].dt.strftime("%Y-%m")))
        
    elif c.REF == 'hadcrut':
        ref = ncr.open_hadcrut(c.DIRECTORY_OBS, c.FILE_OBS, c.VAR)
        ref = ref.sel(time=slice(ctl.time[0].dt.strftime("%Y-%m"),
                                           ref.time[-1].dt.strftime("%Y-%m")))
        
    """
get variables
    """
    ctl = ctl.get(c.VAR)
    sens = sens.get(c.VAR)
    ref = ref.get(c.VAR)
    
    """
season selection, Clim is the yearly mean
    """
    ctl = ncr.dataset_season(ctl, season)
    sens = ncr.dataset_season(sens, season)
    ref = ncr.dataset_season(ref, season)
#%%
    """
bias correction
    """
    #ctl = ctl-(ctl.mean('member').mean('time')-era5_var.mean('time'))    
    #sens = sens-(sens.mean('member').mean('time')-era5_var.mean('time'))
    
    return ctl, sens, ref