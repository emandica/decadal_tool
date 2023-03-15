#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Feb 13 15:42:50 2023

@author: ema
"""
import xarray as xr

import constants as c
import open_dataset as od
import level_selection as ls
import run_bootstrap as rb
import correlation as corr
import data_agg_corr as da_cor

def run_corr(lead_exp, season):
#%%data aggregation
    print('starting data aggregation')

    da_cor.aggr_datasets(lead_exp, season)
    print('data aggregation: OK')
        
    print('starting correlation')
#%%level selection
    ls.level_sel(lead_exp, season)
        
    ctl = xr.open_dataset(c.RUN_DIR+c.NAME_CTL+'_'+c.VAR+'_lead'+str(lead_exp)+'_'+season+'_s'+str(c.T_START)+'_level_'+str(c.PLEV)+'.nc', chunks={'lon':'auto','lat':'auto'})
    sens = xr.open_dataset(c.RUN_DIR+c.NAME_SENS+'_'+c.VAR+'_lead'+str(lead_exp)+'_'+season+'_s'+str(c.T_START)+'_level_'+str(c.PLEV)+'.nc', chunks={'lon':'auto','lat':'auto'})
#%%bootstrap
    rb.bootstrap(ctl,sens,lead_exp,season)
        
    corr.acc(lead_exp, season)
    print('ensemble mean ACC: OK')
        
    corr.acc_significance(lead_exp,season)
    print('correlation significance: OK')
        
    corr.correlation_plot(lead_exp,season)
    print('correlation plot: OK')
    print('correlation: OK')  