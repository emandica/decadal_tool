#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Feb  3 10:02:50 2023

@author: ema
"""
import xarray as xr

import constants as c

def quantile_calculation(lead_exp,season):
    ctl = xr.open_dataset(c.RUN_DIR+c.NAME_CTL+'_'+c.VAR+'_lead'+str(lead_exp)+'_'+season+'_s'+str(c.T_START)+'_bootstrap_standard.nc', chunks={'plev':'auto','lon':'auto','lat':'auto'})
    sens = xr.open_dataset(c.RUN_DIR+c.NAME_SENS+'_'+c.VAR+'_lead'+str(lead_exp)+'_'+season+'_s'+str(c.T_START)+'_bootstrap_standard.nc', chunks={'plev':'auto','lon':'auto','lat':'auto'})
#%%
    """
    bias difference
    """
    delta = sens - ctl
#%%
    """
    quantile calculation
    """
    sign_delta = delta.quantile([0.025,0.05,0.10,0.90,0.95,0.975], dim='bootstrap')
    return sign_delta
