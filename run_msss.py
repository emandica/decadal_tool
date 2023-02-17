#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Feb 13 17:09:05 2023

@author: ema
"""
import xarray as xr

import msss as ms

import data_agg as da
import level_selection as ls

import constants as c

import run_bootstrap as rb


def run_msss(lead_exp, season):
#%%data aggregation
    print('starting data aggregation')
        
    da.aggr_datasets(lead_exp, season)
    print('data aggregation: OK')
        
    print('starting MSSS')
#%%level selection
    ls.level_sel(lead_exp, season)
        
    ctl = xr.open_dataset(c.RUN_DIR+c.NAME_CTL+'_'+c.VAR+'_lead'+str(lead_exp)+'_'+season+'_s'+str(c.T_START)+'_level_'+str(c.PLEV)+'.nc', chunks={'lon':'auto','lat':'auto'})
    sens = xr.open_dataset(c.RUN_DIR+c.NAME_SENS+'_'+c.VAR+'_lead'+str(lead_exp)+'_'+season+'_s'+str(c.T_START)+'_level_'+str(c.PLEV)+'.nc', chunks={'lon':'auto','lat':'auto'})
        
#%%bootstrap
    rb.bootstrap(ctl,sens,lead_exp,season)
        
    ctl.close()
    sens.close()
        
    ms.msss_calc(lead_exp,season)
    print('ensemble mean MSSS: OK')
        
    ms.msss_significance(lead_exp, season)
    print('msss significance: OK')
        
    ms.msss_plot(lead_exp,season)
    print('msss plot: OK')
    print('msss: OK')    