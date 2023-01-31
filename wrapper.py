#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jan 23 16:58:41 2023

@author: ema
"""
import argparse

import generic_routines as cr

import constants as c

import data_agg as da

import run_bootstrap as rb

import mean_bias as mb

import correlation as corr

import plot as pl

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

    season = 'JJA' #args.SEASON
    lead_exp = [3,4] #args.LEAD_LIST

#%%data aggregation    
    print('starting data aggregation')
    
    da.aggr_datasets(lead_exp, season)
    print('data aggregation: OK')
    
#%%bootstrap
    rb.bootstrap(lead_exp,season)
    print('bootstrap: OK')

#%%Mean bias
    if c.M_BIAS:
        print('starting mean bias')        
        
        mb.m_bias(lead_exp, season)
        print('ensemble mean bias: OK')
        
        mb.bias_significance(lead_exp, season)
        print('bias significance: OK')
        
        mb.bias_plot(lead_exp, season)
        print('bias plot: OK')
        print('bias: OK')

                
#%%Correlation
    if c.CORRELATION:    
        print('starting correlation')
        
        corr.acc(lead_exp, season)
        print('ensemble mean ACC: OK')
        
        corr.acc_significance(lead_exp,season)
        print('correlation significance: OK')
        
        corr.correlation_plot(lead_exp,season)
        print('correlation plot: OK')
        print('correlation: OK')    

#%%MSSS        
    #if c.MSSS:
