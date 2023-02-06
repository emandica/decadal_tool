#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Feb  2 11:57:59 2023

@author: ema
"""
import xarray as xr

import numpy as np

import constants as c

import generic_routines as cr

def xsection(lead_exp,season):
    """
    open dataset
    """
    ctl=xr.open_dataset(c.RUN_DIR+c.NAME_CTL+'_'+c.VAR+'_lead'+str(lead_exp)+'_'+season+'_s'+str(c.T_START)+'.nc')
    sens=xr.open_dataset(c.RUN_DIR+c.NAME_SENS+'_'+c.VAR+'_lead'+str(lead_exp)+'_'+season+'_s'+str(c.T_START)+'.nc')
    ref=xr.open_dataset(c.RUN_DIR+c.REF+'_'+c.VAR+'_lead'+str(lead_exp)+'_'+season+'_s'+str(c.T_START)+'.nc')
    
    ctl = cr.cross_section(ctl)
    sens = cr.cross_section(sens)
    ref = cr.cross_section(ref)
    
    """
    time average
    """
    ctl = ctl.mean('time')
    sens = sens.mean('time')
    ref = ref.mean('time')
    """
    save xsection
    """
    ctl.to_netcdf(c.RUN_DIR+c.NAME_CTL+'_'+c.VAR+'_lead'+str(lead_exp)+'_'+season+'_s'+str(c.T_START)+'_xsection.nc')
    sens.to_netcdf(c.RUN_DIR+c.NAME_SENS+'_'+c.VAR+'_lead'+str(lead_exp)+'_'+season+'_s'+str(c.T_START)+'_xsection.nc')    
    ref.to_netcdf(c.RUN_DIR+c.REF+'_'+c.VAR+'_lead'+str(lead_exp)+'_'+season+'_s'+str(c.T_START)+'_xsection.nc')    
    