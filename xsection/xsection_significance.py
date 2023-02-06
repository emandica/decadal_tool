#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Feb  2 12:40:46 2023

@author: ema
"""

import numpy as np

import xarray as xr
from dask.diagnostics import ProgressBar

from scipy import stats

import constants as c

import generic_routines as cr

def xsection_significance(lead_exp,season):
    ctl = xr.open_dataset(c.RUN_DIR+c.NAME_CTL+'_'+c.VAR+'_lead'+str(lead_exp)+'_'+season+'_s'+str(c.T_START)+'.nc', chunks={'plev':'auto','lon':'auto','lat':'auto'}).load()
    sens = xr.open_dataset(c.RUN_DIR+c.NAME_SENS+'_'+c.VAR+'_lead'+str(lead_exp)+'_'+season+'_s'+str(c.T_START)+'.nc', chunks={'plev':'auto','lon':'auto','lat':'auto'}).load()
    ref = xr.open_dataset(c.RUN_DIR+c.REF+'_'+c.VAR+'_lead'+str(lead_exp)+'_'+season+'_s'+str(c.T_START)+'.nc', chunks={'plev':'auto','lon':'auto','lat':'auto'})

    """
    xsection calculation
    """
    ctl = cr.cross_section(ctl)
    sens = cr.cross_section(sens)
    ref = cr.cross_section(ref)    
    
    ctl = ctl.mean('member')
    sens = sens.mean('member')
    
    """
    p value
    """    
    _, p_ctl = stats.ttest_ind(ctl[c.VAR].to_numpy(), ref[c.VAR].to_numpy(), axis=0,equal_var=False)
    _, p_sens = stats.ttest_ind(sens[c.VAR].to_numpy(), ref[c.VAR].to_numpy(), axis=0,equal_var=False)
    
    p_ctl = xr.DataArray(p_ctl, coords=ctl[c.VAR].mean('time').coords)
    
    p_sens = xr.DataArray(p_sens, coords=sens[c.VAR].mean('time').coords)

#%%
    """
    save files
    """
    p_ctl.to_netcdf(c.RUN_DIR+c.NAME_CTL+'_'+c.VAR+'_lead'+str(lead_exp)+'_'+season+'_s'+str(c.T_START)+'_xsection_p_'+c.REF+'.nc')
    p_sens.to_netcdf(c.RUN_DIR+c.NAME_SENS+'_'+c.VAR+'_lead'+str(lead_exp)+'_'+season+'_s'+str(c.T_START)+'_xsection_p_'+c.REF+'.nc')
