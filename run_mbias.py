#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Feb 13 17:03:26 2023

@author: ema
"""
import xarray as xr

import mean_bias as mb

import data_agg as da
import level_selection as ls
import constants as c
import run_bootstrap as rb

def run_mbias(lead_exp, season):
    """
    calculate variable mean bias
    """
#%%data aggregation
    print('starting data aggregation')

    da.aggr_datasets(lead_exp, season)
    print('data aggregation: OK')

    print('starting mean bias')
    ls.level_sel(lead_exp, season)

#%%bootstrap
    ctl_path = c.RUN_DIR+c.NAME_CTL+'_'+c.VAR+'_lead'+str(lead_exp)+'_'+season+'_s'+str(c.T_START)+'_level_'+str(c.PLEV)+'.nc'
    ctl = xr.open_dataset(ctl_path, chunks={'lon':'auto','lat':'auto'})

    sens_path = c.RUN_DIR+c.NAME_SENS+'_'+c.VAR+'_lead'+str(lead_exp)+'_'+season+'_s'+str(c.T_START)+'_level_'+str(c.PLEV)+'.nc'
    sens = xr.open_dataset(sens_path, chunks={'lon':'auto','lat':'auto'})

    ctl = ctl.mean('time')
    sens = sens.mean('time')

    rb.bootstrap(ctl,sens,lead_exp,season)

    ctl.close()
    sens.close()
    print('bootstrap: OK')

    mb.bias_significance(lead_exp, season)
    print('bias significance: OK')

    mb.bias_plot(lead_exp, season)
    print('bias plot: OK')
    print('bias: OK')
