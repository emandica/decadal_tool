#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Feb 13 16:58:51 2023

@author: ema
"""

import xarray as xr

import xsection as xs

import data_agg as da
import constants as c
import run_bootstrap as rb

def run_xsec(lead_exp, season):
#%%data aggregation 
    print('starting data aggregation')

    da.aggr_datasets(lead_exp, season)
    print('data aggregation: OK')
        
        
    print('starting crossection')
    xs.xsection(lead_exp, season)
    print('ensemble mean xsection: OK')
        
#%%bootstrap
    ctl = xr.open_dataset(c.RUN_DIR+c.NAME_CTL+'_'+c.VAR+'_lead'+str(lead_exp)+'_'+season+'_s'+str(c.T_START)+'_xsection.nc', chunks={'lon':'auto','lat':'auto'})
    sens = xr.open_dataset(c.RUN_DIR+c.NAME_SENS+'_'+c.VAR+'_lead'+str(lead_exp)+'_'+season+'_s'+str(c.T_START)+'_xsection.nc', chunks={'lon':'auto','lat':'auto'})
            
    rb.bootstrap(ctl,sens,lead_exp,season)
        
    ctl.close()
    sens.close()
    print('bootstrap: OK')
        
    xs.xsection_significance(lead_exp, season)
    print('xsection significance: OK')

    xs.xsection_plot(lead_exp, season)      
    print('xsection plt: OK')
    print('xsection: OK')