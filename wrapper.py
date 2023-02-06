#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jan 23 16:58:41 2023

@author: ema
"""
import argparse
import xarray as xr

import mean_bias as mb
import correlation as corr
import xsection as xs

import constants as c

import data_agg as da

import run_bootstrap as rb

import level_selection as ls
#%%
# entry point for the program
if __name__ == '__main__':
###############################################################################
###constants definitiontion
    CLI=argparse.ArgumentParser()
    CLI.add_argument(
        "--LEAD_LIST",
        nargs="*",
        type=int,
        default=[],
        )
    
    CLI.add_argument(
            "-SEASON",
            #nargs="*",
            type=str,
            default=[],
            )

    args=CLI.parse_args()

    season = 'DJF' #args.SEASON
    lead_exp = [3,4] #args.LEAD_LIST

#%%data aggregation 
    print('starting data aggregation')
    
    da.aggr_datasets(lead_exp, season)
    print('data aggregation: OK')

#%%xsection
    if c.XSECT:
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
        
#%%Mean bias
    if c.M_BIAS:
        
        print('starting mean bias')        
        ls.level_sel(lead_exp, season)
        #mb.m_bias(lead_exp, season)
        #print('ensemble mean bias: OK')
       
        #%%bootstrap
        ctl = xr.open_dataset(c.RUN_DIR+c.NAME_CTL+'_'+c.VAR+'_lead'+str(lead_exp)+'_'+season+'_s'+str(c.T_START)+'_level_'+str(c.PLEV)+'.nc', chunks={'lon':'auto','lat':'auto'})
        sens = xr.open_dataset(c.RUN_DIR+c.NAME_SENS+'_'+c.VAR+'_lead'+str(lead_exp)+'_'+season+'_s'+str(c.T_START)+'_level_'+str(c.PLEV)+'.nc', chunks={'lon':'auto','lat':'auto'})
        
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

                
#%%Correlation
    if c.CORRELATION:
        
        print('starting correlation')

        ls.level_sel(lead_exp, season)
        
        #%%bootstrap
        ctl = xr.open_dataset(c.RUN_DIR+c.NAME_CTL+'_'+c.VAR+'_lead'+str(lead_exp)+'_'+season+'_s'+str(c.T_START)+'_level_'+str(c.PLEV)+'.nc', chunks={'lon':'auto','lat':'auto'})
        sens = xr.open_dataset(c.RUN_DIR+c.NAME_SENS+'_'+c.VAR+'_lead'+str(lead_exp)+'_'+season+'_s'+str(c.T_START)+'_level_'+str(c.PLEV)+'.nc', chunks={'lon':'auto','lat':'auto'})
        
        rb.bootstrap(ctl,sens,lead_exp,season)
        
        ctl.close()
        sens.close()
        
        corr.acc(lead_exp, season)
        print('ensemble mean ACC: OK')
        
        corr.acc_significance(lead_exp,season)
        print('correlation significance: OK')
        
        corr.correlation_plot(lead_exp,season)
        print('correlation plot: OK')
        print('correlation: OK')    

#%%MSSS        
    #if c.MSSS:
        