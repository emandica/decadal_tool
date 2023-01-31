#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jan 27 13:15:40 2023

@author: ema
"""
import xarray as xr

import constants as c

#import confess_routines as cr

def m_bias(lead_exp,season):
    """
    open dataset
    """
    ctl=xr.open_dataset(c.RUN_DIR+c.NAME_CTL+'_'+c.VAR+'_lead'+str(lead_exp)+'_'+season+'_s'+str(c.T_START)+'.nc')
    sens=xr.open_dataset(c.RUN_DIR+c.NAME_SENS+'_'+c.VAR+'_lead'+str(lead_exp)+'_'+season+'_s'+str(c.T_START)+'.nc')
    ref=xr.open_dataset(c.RUN_DIR+c.REF+'_'+c.VAR+'_lead'+str(lead_exp)+'_'+season+'_s'+str(c.T_START)+'.nc')
    
    """
    ensemble mean
    """
    ctl = ctl.mean('member')
    sens = sens.mean('member')
    """
    time mean
    """
    ctl = ctl.mean('time')
    sens = sens.mean('time')
    ref = ref.mean('time')
    
    ctl = ctl - ref
    sens = sens - ref
    
    """
    save mean bias
    """
    ctl.to_netcdf(c.RUN_DIR+c.NAME_CTL+'_'+c.VAR+'_lead'+str(lead_exp)+'_'+season+'_s'+str(c.T_START)+'_MBIAS.nc')
    sens.to_netcdf(c.RUN_DIR+c.NAME_SENS+'_'+c.VAR+'_lead'+str(lead_exp)+'_'+season+'_s'+str(c.T_START)+'_MBIAS.nc')