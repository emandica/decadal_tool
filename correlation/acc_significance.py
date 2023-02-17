#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Jan 22 17:30:29 2023

@author: ema
"""

import xarray as xr
import xskillscore as xs

import constants as c

###############################################################################
def acc_significance(lead_exp,season):
    
    if c.BOOTSTRAP == 'block':
        ctl = xr.open_dataset(c.RUN_DIR+c.NAME_CTL+'_'+c.VAR+'_lead'+str(lead_exp)+'_'+season+'_s'+str(c.T_START)+'_level_'+str(c.PLEV)+'.nc', chunks={'lon':'auto','lat':'auto'})
        sens = xr.open_dataset(c.RUN_DIR+c.NAME_SENS+'_'+c.VAR+'_lead'+str(lead_exp)+'_'+season+'_s'+str(c.T_START)+'_level_'+str(c.PLEV)+'.nc', chunks={'lon':'auto','lat':'auto'})
        ref = xr.open_dataset(c.RUN_DIR+c.REF+'_'+c.VAR+'_lead'+str(lead_exp)+'_'+season+'_s'+str(c.T_START)+'_level_'+str(c.PLEV)+'.nc', chunks={'lon':'auto','lat':'auto'})

#%%    
        """
        Anomaly Correlations coefficient
        """
        corr_ctl = xr.corr(ctl[c.VAR], ref[c.VAR], dim='time')
        corr_sens = xr.corr(sens[c.VAR], ref[c.VAR], dim='time')

#%%
        """
        save files
        """
        corr_ctl.to_netcdf(c.RUN_DIR+c.NAME_CTL+'_'+c.VAR+'_lead'+str(lead_exp)+'_'+season+'_s'+str(c.T_START)+'_bootstrap_correlation_'+c.REF+'.nc')
        corr_sens.to_netcdf(c.RUN_DIR+c.NAME_SENS+'_'+c.VAR+'_lead'+str(lead_exp)+'_'+season+'_s'+str(c.T_START)+'_bootstrap_correlation_'+c.REF+'.nc')


    elif c.BOOTSTRAP == 'standard':
        
        ctl = xr.open_dataset(c.RUN_DIR+c.NAME_CTL+'_'+c.VAR+'_lead'+str(lead_exp)+'_'+season+'_s'+str(c.T_START)+'_level_'+str(c.PLEV)+'.nc', chunks={'lon':'auto','lat':'auto'})
        sens = xr.open_dataset(c.RUN_DIR+c.NAME_SENS+'_'+c.VAR+'_lead'+str(lead_exp)+'_'+season+'_s'+str(c.T_START)+'_level_'+str(c.PLEV)+'.nc', chunks={'lon':'auto','lat':'auto'})
        delta1 = xr.open_dataset(c.RUN_DIR+c.NAME_CTL+'_'+c.VAR+'_lead'+str(lead_exp)+'_'+season+'_s'+str(c.T_START)+'_bootstrap_standard.nc', chunks={'lon':'auto','lat':'auto'})
        delta2 = xr.open_dataset(c.RUN_DIR+c.NAME_SENS+'_'+c.VAR+'_lead'+str(lead_exp)+'_'+season+'_s'+str(c.T_START)+'_bootstrap_standard.nc', chunks={'lon':'auto','lat':'auto'})
        ref = xr.open_dataset(c.RUN_DIR+c.REF+'_'+c.VAR+'_lead'+str(lead_exp)+'_'+season+'_s'+str(c.T_START)+'_level_'+str(c.PLEV)+'.nc', chunks={'lon':'auto','lat':'auto'})
        
        """
        delta ACC
        """
        delta1_corr = xr.corr(delta1[c.VAR],ref[c.VAR],dim='time')
        delta2_corr = xr.corr(delta2[c.VAR],ref[c.VAR],dim='time')
        
        """
        Correlation p value
        """
        ctl =ctl.mean('member')
        sens =sens.mean('member')
        
        corr_ctl = xs.pearson_r_eff_p_value(ctl, ref, dim='time',keep_attrs=True)
        corr_sens = xs.pearson_r_eff_p_value(sens, ref, dim='time',keep_attrs=True)
        
        """
        save file
        """
        corr_ctl.to_netcdf(c.RUN_DIR+c.NAME_CTL+'_'+c.VAR+'_lead'+str(lead_exp)+'_'+season+'_s'+str(c.T_START)+'_correlation_p_'+c.REF+'.nc')
        corr_sens.to_netcdf(c.RUN_DIR+c.NAME_SENS+'_'+c.VAR+'_lead'+str(lead_exp)+'_'+season+'_s'+str(c.T_START)+'_correlation_p_'+c.REF+'.nc')
        delta1_corr.to_netcdf(c.RUN_DIR+c.NAME_CTL+'_'+c.VAR+'_lead'+str(lead_exp)+'_'+season+'_s'+str(c.T_START)+'_bootstrap_standard_correlation_'+c.REF+'.nc')
        delta2_corr.to_netcdf(c.RUN_DIR+c.NAME_SENS+'_'+c.VAR+'_lead'+str(lead_exp)+'_'+season+'_s'+str(c.T_START)+'_bootstrap_standard_correlation_'+c.REF+'.nc')
###############################################################################    
#%%
#    """
#    land-sea mask
#    """
#    land_ctl, oce_ctl = ncr.land_sea_mask(corr_ctl)
#    land_sens, oce_sens = ncr.land_sea_mask(corr_sens)