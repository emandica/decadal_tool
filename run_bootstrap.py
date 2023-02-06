#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Jan 22 16:52:29 2023

@author: ema
"""

import xarray as xr

import constants as c
import b_bootstrap as b
import generic_routines as cr
###############################################################################
def bootstrap(ctl, sens, lead_exp, season):
    
    if c.BOOTSTRAP=='block':
        ctl = b.block_bootstrap(ctl,1000)
        sens = b.block_bootstrap(sens,1000)
    
        ctl.to_netcdf(c.RUN_DIR+c.NAME_CTL+'_'+c.VAR+'_lead'+str(lead_exp)+'_'+season+'_s'+str(c.T_START)+'_bootstrap.nc')
        ctl.close()
        
        sens.to_netcdf(c.RUN_DIR+c.NAME_SENS+'_'+c.VAR+'_lead'+str(lead_exp)+'_'+season+'_s'+str(c.T_START)+'_bootstrap.nc')
        sens.close()
        
    elif c.BOOTSTRAP=='standard':        
        
        cr.exp_concat(ctl, sens, 'member', lead_exp, season)
        
        dset= xr.open_dataset(c.RUN_DIR+c.NAME_CTL+'_'+c.VAR+'_lead'+str(lead_exp)+'_'+season+'_s'+str(c.T_START)+'_concat.nc', chunks={'lon':'auto', 'lat':'auto'})
        
        dset1= b.standard_bootstrap(dset,1000)
        dset1.to_netcdf(c.RUN_DIR+c.NAME_CTL+'_'+c.VAR+'_lead'+str(lead_exp)+'_'+season+'_s'+str(c.T_START)+'_bootstrap_standard.nc')
        dset1.close()
        
        dset2= b.standard_bootstrap(dset,1000)
        dset2.to_netcdf(c.RUN_DIR+c.NAME_SENS+'_'+c.VAR+'_lead'+str(lead_exp)+'_'+season+'_s'+str(c.T_START)+'_bootstrap_standard.nc')
        dset2.close()
    