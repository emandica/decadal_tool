#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Mar 22 10:05:51 2023

@author: ema
"""

import numpy as np

from eofs.xarray import Eof

import xarray as xr

import os.path

import constants as c
import data_agg as da
import level_selection as ls
import run_bootstrap as rb

def pcs_check(var_eofs,var_pcs):
    
    if var_eofs[0].sel(**{'lat':slice(63,70),'lon':slice(-25,-16)}).mean('lat').mean('lon') > 0:
        var_eofs=-var_eofs
        var_pcs=-var_pcs
        
    return var_eofs, var_pcs
        
    
def nao_index_eofs(var):
    var.coords['lon']=(var.coords['lon'] + 180) % 360 - 180
    var=var.sortby(var.lon)
    LatIndexer, LonIndexer = 'lat', 'lon'
    var = var.sel(**{LatIndexer: slice(20,80),
                     LonIndexer: slice(-90,40)})
    
    wgts = np.sqrt(np.cos(np.deg2rad(var.lat.values)))[...,np.newaxis]
    solver = Eof(var, weights=wgts)
    pcs = solver.pcs(npcs=1, pcscaling=1)
    eofs = solver.eofsAsCorrelation(neofs=1)
    
    return pcs, eofs

def nao_loop(dset):
    nao_pcs=[]
    nao_eofs=[]
    for i in range(1000):
        dset_mem = dset.isel({'bootstrap':i})
        anci1,anci2=nao_index_eofs(dset_mem)
        nao_pcs.append(anci1)
        nao_eofs.append(anci2)
    nao_pcs = xr.concat(nao_pcs, dim='bootstrap')
    nao_eofs = xr.concat(nao_eofs, dim='bootstrap')
    nao_pcs = nao_pcs.squeeze()
    nao_eofs = nao_eofs.squeeze()
    return nao_pcs, nao_eofs

def NAO_bootstrap(lead_exp,season):
    if os.path.isfile(c.RUN_DIR+c.NAME_CTL+'_'+c.VAR+'_lead'+str(lead_exp)+'_'+season+'_s'+str(c.T_START)+'.nc') == False:
        da.aggr_datasets(lead_exp, season)
    if os.path.isfile(c.RUN_DIR+c.NAME_CTL+'_'+c.VAR+'_lead'+str(lead_exp)+'_'+season+'_s'+str(c.T_START)+'_level_'+str(c.PLEV)+'.nc') == False:            
        ls.level_sel(lead_exp, season)
    if  os.path.isfile(c.RUN_DIR+c.NAME_CTL+'_'+c.VAR+'_lead'+str(lead_exp)+'_'+season+'_s'+str(c.T_START)+'bootstrap_standard.nc') == False:       
        ctl = xr.open_dataset(c.RUN_DIR+c.NAME_CTL+'_'+c.VAR+'_lead'+str(lead_exp)+'_'+season+'_s'+str(c.T_START)+'_level_'+str(c.PLEV)+'.nc', chunks={'lon':'auto','lat':'auto'})
        sens = xr.open_dataset(c.RUN_DIR+c.NAME_SENS+'_'+c.VAR+'_lead'+str(lead_exp)+'_'+season+'_s'+str(c.T_START)+'_level_'+str(c.PLEV)+'.nc', chunks={'lon':'auto','lat':'auto'})
        rb.bootstrap(ctl,sens,lead_exp,season)
        
    delta1 = xr.open_dataset(c.RUN_DIR+c.NAME_CTL+'_'+c.VAR+'_lead'+str(lead_exp)+'_'+season+'_s'+str(c.T_START)+'_bootstrap_standard.nc')
    delta2 = xr.open_dataset(c.RUN_DIR+c.NAME_SENS+'_'+c.VAR+'_lead'+str(lead_exp)+'_'+season+'_s'+str(c.T_START)+'_bootstrap_standard.nc')        
    era_pcs = xr.open_dataset(c.RUN_DIR+c.REF+'_nao_index_lead'+str(lead_exp)+'_'+season+'.nc')
    
    nao_d1_pcs, nao_d1_eofs=nao_loop(delta1[c.VAR])
    nao_d2_pcs, nao_d2_eofs=nao_loop(delta2[c.VAR])

    d1_pcs_cor = xr.corr(nao_d1_pcs,era_pcs.pcs,dim='time')
    d2_pcs_cor = xr.corr(nao_d2_pcs,era_pcs.pcs,dim='time')

    delta = d1_pcs_cor -d2_pcs_cor
    
    quant = delta.quantile([0.025,0.05,0.10,0.90,0.95,0.975], dim='bootstrap')
