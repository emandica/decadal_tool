#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jan 18 11:47:05 2023

@author: ema
"""

import generic_routines as ncr
import constants as c

def open_dataset_and_seasonal_selection(lead_exp, season):
    
    """
Load datasets
    """    
    ctl = ncr.lead_aggregation(lead_exp, c.VAR, c.CTL_NUMBER, c.NAME_CTL,
                               c.T_START,
                               no_member=c.no_member,
                               )
    print('aggregation ctl: OK')
    
    sens = ncr.lead_aggregation(lead_exp, c.VAR, c.SENS_NUMBER, c.NAME_SENS,
                                c.T_START,
                                no_member=c.no_member,
                                )
    print('aggregation sens: OK')

    if c.REF=='era':
        ref = ncr.open_era5(c.DIRECTORY_ERA, c.FILE_ERA, c.VAR)
        ref = ref.sel(time=slice(ctl.time[0],
                                 ref.time[-1]))
        
    elif c.REF == 'hadcrut':
        ref = ncr.open_hadcrut(c.DIRECTORY_OBS, c.FILE_OBS, c.VAR)
        ref = ref.sel(time=slice(ctl.time[0].dt.strftime("%Y-%m"),
                                 ref.time[-1].dt.strftime("%Y-%m")))
    print('reference: OK')        
        
    """
season selection, Clim is the yearly mean
    """
    ctl_s = ncr.dataset_season(ctl[c.VAR], season)
    sens_s = ncr.dataset_season(sens[c.VAR], season)
    ref_s = ncr.dataset_season(ref[c.VAR], season)
    
    return ctl_s, sens_s, ref_s