#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jan 27 14:04:55 2023

@author: ema
"""
import numpy as np

import xarray as xr
#import xarray.apply_ufuncs as xrf

from scipy import stats

import constants as c

###############################################################################
def bias_significance(lead_exp,season):
    ctl = xr.open_dataset(c.RUN_DIR+c.NAME_CTL+'_'+c.VAR+'_lead'+str(lead_exp)+'_'+season+'_s'+str(c.T_START)+'.nc', chunks={'lon':'auto','lat':'auto'}).load()
    sens = xr.open_dataset(c.RUN_DIR+c.NAME_SENS+'_'+c.VAR+'_lead'+str(lead_exp)+'_'+season+'_s'+str(c.T_START)+'.nc', chunks={'lon':'auto','lat':'auto'}).load()
    delta1 = xr.open_dataset(c.RUN_DIR+c.NAME_CTL+'_'+c.VAR+'_lead'+str(lead_exp)+'_'+season+'_s'+str(c.T_START)+'_bootstrap_standard.nc', chunks={'lon':'auto','lat':'auto'})
    delta2 = xr.open_dataset(c.RUN_DIR+c.NAME_SENS+'_'+c.VAR+'_lead'+str(lead_exp)+'_'+season+'_s'+str(c.T_START)+'_bootstrap_standard.nc', chunks={'lon':'auto','lat':'auto'})
    ref = xr.open_dataset(c.RUN_DIR+c.REF+'_'+c.VAR+'_lead'+str(lead_exp)+'_'+season+'_s'+str(c.T_START)+'.nc', chunks={'lon':'auto','lat':'auto'})
    
    """
    mean bias
    """
    mb_ctl = delta1.mean('time') - ref.mean('time')
    mb_sens = delta2.mean('time') - ref.mean('time')
    
    """
    ensemble mean
    """
    ctl=ctl.mean('member')
    sens=sens.mean('member')
    
    """
    Bias p value
    """
    
    _, p_ctl = stats.ttest_ind(ctl[c.VAR].to_numpy(), ref[c.VAR].to_numpy(), axis=0,equal_var=False)
    _, p_sens = stats.ttest_ind(sens[c.VAR].to_numpy(), ref[c.VAR].to_numpy(), axis=0,equal_var=False)
    
    #t_ctl = (ctl.mean('time')-ref.mean('time')/xr.apply_ufunc(np.sqrt,(ctl.var(dim='time',ddof=1))/ float(len(ctl['time']))))
    #t_sens = (sens.mean('time')-ref.mean('time')/xr.apply_ufunc(np.sqrt,(sens.var(dim='time',ddof=1))/ float(len(sens['time']))))
    
   # prob_ctl = stats.distributions.t.sf(xr.apply_ufunc(np.fabs,t_ctl[c.VAR].load()), len(ctl['time'])-1) * 2
    p_ctl = xr.DataArray(p_ctl, coords=ctl[c.VAR].mean('time').coords)
    
    #prob_sens = stats.distributions.t.sf(xr.apply_ufunc(np.fabs,t_sens[c.VAR].load()), len(sens['time'])-1) * 2
    p_sens = xr.DataArray(p_sens, coords=sens[c.VAR].mean('time').coords)
#%%
    """
    save files
    """
    mb_ctl.to_netcdf(c.RUN_DIR+c.NAME_CTL+'_'+c.VAR+'_lead'+str(lead_exp)+'_'+season+'_s'+str(c.T_START)+'_bootstrap_mbias_'+c.REF+'.nc')
    mb_sens.to_netcdf(c.RUN_DIR+c.NAME_SENS+'_'+c.VAR+'_lead'+str(lead_exp)+'_'+season+'_s'+str(c.T_START)+'_bootstrap_mbias_'+c.REF+'.nc')
    p_ctl.to_netcdf(c.RUN_DIR+c.NAME_CTL+'_'+c.VAR+'_lead'+str(lead_exp)+'_'+season+'_s'+str(c.T_START)+'_bias_p_'+c.REF+'.nc')
    p_sens.to_netcdf(c.RUN_DIR+c.NAME_SENS+'_'+c.VAR+'_lead'+str(lead_exp)+'_'+season+'_s'+str(c.T_START)+'_bias_p_'+c.REF+'.nc')