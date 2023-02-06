#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Nov 29 20:43:24 2022

@author: ema
"""
import xarray as xr

import constants as c

###############################################################################
def acc(lead_exp, season):
    """
    open dataset
    """
    ctl=xr.open_dataset(c.RUN_DIR+c.NAME_CTL+'_'+c.VAR+'_lead'+str(lead_exp)+'_'+season+'_s'+str(c.T_START)+'_level_'+str(c.PLEV)+'.nc')
    sens=xr.open_dataset(c.RUN_DIR+c.NAME_SENS+'_'+c.VAR+'_lead'+str(lead_exp)+'_'+season+'_s'+str(c.T_START)+'_level_'+str(c.PLEV)+'.nc')
    ref=xr.open_dataset(c.RUN_DIR+c.REF+'_'+c.VAR+'_lead'+str(lead_exp)+'_'+season+'_s'+str(c.T_START)+'_level_'+str(c.PLEV)+'.nc')
    
#%%
    """
    ACC
    """
    ctl =ctl.mean('member')
    sens =sens.mean('member')
    
    ctl = xr.corr(ctl[c.VAR], ref[c.VAR], dim='time')    
    sens = xr.corr(sens[c.VAR], ref[c.VAR], dim='time')

    """
    save ACC
    """
    ctl.to_netcdf(c.RUN_DIR+c.NAME_CTL+'_'+c.VAR+'_lead'+str(lead_exp)+'_'+season+'_s'+str(c.T_START)+'_ACC.nc')
    sens.to_netcdf(c.RUN_DIR+c.NAME_SENS+'_'+c.VAR+'_lead'+str(lead_exp)+'_'+season+'_s'+str(c.T_START)+'_ACC.nc')
