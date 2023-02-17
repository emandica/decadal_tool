#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Feb  9 12:34:20 2023

@author: ema
"""

import xarray as xr

from scipy import stats

import constants as c

import generic_routines as gr

def msss_significance(lead_exp,season):
    
    if c.BOOTSTRAP == 'standard':
            
        ctl = xr.open_dataset(c.RUN_DIR+c.NAME_CTL+'_'+c.VAR+'_lead'+str(lead_exp)+'_'+season+'_s'+str(c.T_START)+'_level_'+str(c.PLEV)+'.nc', chunks={'lon':'auto','lat':'auto'})
        sens = xr.open_dataset(c.RUN_DIR+c.NAME_SENS+'_'+c.VAR+'_lead'+str(lead_exp)+'_'+season+'_s'+str(c.T_START)+'_level_'+str(c.PLEV)+'.nc', chunks={'lon':'auto','lat':'auto'})
        delta1 = xr.open_dataset(c.RUN_DIR+c.NAME_CTL+'_'+c.VAR+'_lead'+str(lead_exp)+'_'+season+'_s'+str(c.T_START)+'_bootstrap_standard.nc', chunks={'lon':'auto','lat':'auto'})
        delta2 = xr.open_dataset(c.RUN_DIR+c.NAME_SENS+'_'+c.VAR+'_lead'+str(lead_exp)+'_'+season+'_s'+str(c.T_START)+'_bootstrap_standard.nc', chunks={'lon':'auto','lat':'auto'})
        ref = xr.open_dataset(c.RUN_DIR+c.REF+'_'+c.VAR+'_lead'+str(lead_exp)+'_'+season+'_s'+str(c.T_START)+'_level_'+str(c.PLEV)+'.nc', chunks={'lon':'auto','lat':'auto'})

        """
        mse
        """
        delta1 = gr.mse(delta1,ref)
        delta2 = gr.mse(delta2,ref)

        ctl = gr.mse(ctl, ref)
        sens = gr.mse(sens, ref)
        zero = gr.mse(ref.mean('time'),ref)
        
        """
        p values
        """
        _,p_ctl = stats.ttest_1samp(ctl[c.VAR].to_numpy(), zero[c.VAR].to_numpy())
        p_ctl = xr.DataArray(p_ctl, coords=ctl[c.VAR].mean('member').coords)
        
        _,p_sens = stats.ttest_1samp(sens[c.VAR].to_numpy(), zero[c.VAR].to_numpy())
        p_sens = xr.DataArray(p_sens, coords=sens[c.VAR].mean('member').coords)
        
        #%%
        """
        save files
        """
        delta1.to_netcdf(c.RUN_DIR+c.NAME_CTL+'_'+c.VAR+'_lead'+str(lead_exp)+'_'+season+'_s'+str(c.T_START)+'_bootstrap_mse_'+c.REF+'.nc')
        delta2.to_netcdf(c.RUN_DIR+c.NAME_SENS+'_'+c.VAR+'_lead'+str(lead_exp)+'_'+season+'_s'+str(c.T_START)+'_bootstrap_mse_'+c.REF+'.nc')
        
        ctl.to_netcdf(c.RUN_DIR+c.NAME_CTL+'_'+c.VAR+'_lead'+str(lead_exp)+'_'+season+'_s'+str(c.T_START)+'_mse.nc')
        sens.to_netcdf(c.RUN_DIR+c.NAME_SENS+'_'+c.VAR+'_lead'+str(lead_exp)+'_'+season+'_s'+str(c.T_START)+'_mse.nc')
        zero.to_netcdf(c.RUN_DIR+c.REF+'_'+c.VAR+'_lead'+str(lead_exp)+'_'+season+'_s'+str(c.T_START)+'_mse.nc')
        p_ctl.to_netcdf(c.RUN_DIR+c.NAME_CTL+'_'+c.VAR+'_lead'+str(lead_exp)+'_'+season+'_s'+str(c.T_START)+'_mse_p_'+c.REF+'.nc')
        p_sens.to_netcdf(c.RUN_DIR+c.NAME_SENS+'_'+c.VAR+'_lead'+str(lead_exp)+'_'+season+'_s'+str(c.T_START)+'_mse_p_'+c.REF+'.nc')
            
        
