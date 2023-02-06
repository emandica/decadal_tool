#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jan 27 14:04:55 2023

@author: ema
"""
import numpy as np

import xarray as xr

from scipy import stats

import constants as c

###############################################################################
def bias_significance(lead_exp,season):
    ctl = xr.open_dataset(c.RUN_DIR+c.NAME_CTL+'_'+c.VAR+'_lead'+str(lead_exp)+'_'+season+'_s'+str(c.T_START)+'_level_'+str(c.PLEV)+'.nc', chunks={'lon':'auto','lat':'auto'})
    sens = xr.open_dataset(c.RUN_DIR+c.NAME_SENS+'_'+c.VAR+'_lead'+str(lead_exp)+'_'+season+'_s'+str(c.T_START)+'_level_'+str(c.PLEV)+'.nc', chunks={'lon':'auto','lat':'auto'})
    delta1 = xr.open_dataset(c.RUN_DIR+c.NAME_CTL+'_'+c.VAR+'_lead'+str(lead_exp)+'_'+season+'_s'+str(c.T_START)+'_bootstrap_standard.nc', chunks={'lon':'auto','lat':'auto'})
    delta2 = xr.open_dataset(c.RUN_DIR+c.NAME_SENS+'_'+c.VAR+'_lead'+str(lead_exp)+'_'+season+'_s'+str(c.T_START)+'_bootstrap_standard.nc', chunks={'lon':'auto','lat':'auto'})
    ref = xr.open_dataset(c.RUN_DIR+c.REF+'_'+c.VAR+'_lead'+str(lead_exp)+'_'+season+'_s'+str(c.T_START)+'_level_'+str(c.PLEV)+'.nc', chunks={'lon':'auto','lat':'auto'})

    """
    mean bias
    """
    mb_ctl = delta1 - ref.mean('time')
    mb_sens = delta2 - ref.mean('time')
    
    """
    ensemble mean
    """
    ctl=ctl.mean('member')
    sens=sens.mean('member')
    
    """
    Bias p value
    """    
    _, p_ctl = stats.ttest_ind(ctl[c.VAR].to_numpy(), ref[c.VAR].to_numpy(), axis=0,equal_var=False)    
    p_ctl = xr.DataArray(p_ctl, coords=ctl[c.VAR].mean('time').coords)
    ctl = ctl.mean('time')
    
    _, p_sens = stats.ttest_ind(sens[c.VAR].to_numpy(), ref[c.VAR].to_numpy(), axis=0,equal_var=False)
    p_sens = xr.DataArray(p_sens, coords=sens[c.VAR].mean('time').coords)
    sens = sens.mean('time')
#%%
    """
    save files
    """
    mb_ctl.to_netcdf(c.RUN_DIR+c.NAME_CTL+'_'+c.VAR+'_lead'+str(lead_exp)+'_'+season+'_s'+str(c.T_START)+'_bootstrap_mbias_'+c.REF+'.nc')
    mb_sens.to_netcdf(c.RUN_DIR+c.NAME_SENS+'_'+c.VAR+'_lead'+str(lead_exp)+'_'+season+'_s'+str(c.T_START)+'_bootstrap_mbias_'+c.REF+'.nc')
    
    ctl.to_netcdf(c.RUN_DIR+c.NAME_CTL+'_'+c.VAR+'_lead'+str(lead_exp)+'_'+season+'_s'+str(c.T_START)+'_MBIAS.nc')
    sens.to_netcdf(c.RUN_DIR+c.NAME_SENS+'_'+c.VAR+'_lead'+str(lead_exp)+'_'+season+'_s'+str(c.T_START)+'_MBIAS.nc')
    p_ctl.to_netcdf(c.RUN_DIR+c.NAME_CTL+'_'+c.VAR+'_lead'+str(lead_exp)+'_'+season+'_s'+str(c.T_START)+'_bias_p_'+c.REF+'.nc')
    p_sens.to_netcdf(c.RUN_DIR+c.NAME_SENS+'_'+c.VAR+'_lead'+str(lead_exp)+'_'+season+'_s'+str(c.T_START)+'_bias_p_'+c.REF+'.nc')
