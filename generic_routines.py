#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Jan 21 19:21:11 2023

@author: ema
"""
import numpy as np
import xarray as xr
import pandas as pd

import constants as c

###############################################################################
###FUNCTIONS TO OPEN FILES###
###############################################################################
def open_dataset(path, lead, member_number, time_start, no_member=False):
    """opens netcdf files. Separate ensamble members are concatenate
    Args:
        path (str): path to data (data should be aggregated by lead time)
        lead (int): lead year
        member_number (int): number of ensemble members
        time_start(): First year
        time_end(): Last year
        no_member(int): ignore selected member
    """
    ancillary = []  #variable to collect single members
    for member in range(1, member_number+1):
        if member in no_member:
            print('member '+str(member)+' is ignored')
        else:
            name = "r"+str(member)+"_lead_"+str(lead)+"_r.nc"
            data = xr.open_dataset(path+name,chunks={'time':'auto','lon':'auto','lat':'auto'})
            ancillary.append(data)
            print('member '+str(member))
    data= xr.concat(ancillary, dim='member')
    data = data.sel(time=slice(time_start, data.time[-1]))
    return data

###############################################################################
def lead_aggregation(lead, var, number, name,t_start, no_member):
    """ aggregates the lead times
    Args:
        lead(): lead year
        VAR(): model variable
        NUMBER(): number of ensemble members
        NAME(): exp name
        no_member(): member to ignore
    """"""BOOTSTRAP FOR DECADAL FORECASTS"""
    ancillary=[]
    flag = 0
    for sub_lead in lead:
        
        time_start = t_start+sub_lead
        
        date_start = str(time_start)+"-11"
        
        path = c.FILE_DIR+name+"/monthly/"+var+"/lead_"+str(sub_lead)+"/"

        dset = open_dataset(path, sub_lead, number,
                            time_start=date_start,
                            no_member=no_member)
        
        if flag == 0:
            time_ax = pd.date_range(start=str(date_start), end=pd.to_datetime(dset.time.dt.strftime("%Y-%m"))[-1], freq='MS')
            flag = 1
            
        dset['time']=time_ax
        
        ancillary.append(dset)

    dset = xr.concat(ancillary, dim='lead')
    return dset

###############################################################################
def open_era5(directory, file, variable):
    """converts era variables in CMOR names
    Args:
        diectory():
        file():
        variable():
    """
    era5_var = xr.open_dataset(directory+file,chunks={'time':'auto','lon':'auto','lat':'auto'})
    if variable == 'pr':
        era5_var = era5_var.rename({'tp': 'pr'})
        era5_var = era5_var*1000
    elif variable == 'ua':
        era5_var = era5_var.rename({'u': 'ua'})
    elif variable =='tas':
        era5_var = era5_var.rename({'2t': 'tas'})
    elif variable == 'z':
        era5_var = era5_var.rename({'zg': 'z'})
    elif variable == 'psl':
        era5_var = era5_var.rename({'msl': 'psl'})
    return era5_var

###############################################################################
def open_hadcrut(directory, file, variable):
    """converts era variables in CMOR names
    Args:
        diectory():
        file():
        variable():
    """
    ref_var = xr.open_dataset(directory+file)
    ref_var = ref_var.rename({'tas_mean': 'tas'})
    return ref_var

###############################################################################
def open_dataset_and_seasonal_selection(lead_exp, season):
    
    """
Load datasets
    """    
    ctl = lead_aggregation(lead_exp, c.VAR, c.CTL_NUMBER, c.NAME_CTL,
                               c.T_START,c.T_END,
                               no_member=c.no_member,
                               )

    sens = lead_aggregation(lead_exp, c.VAR, c.SENS_NUMBER, c.NAME_SENS,
                                c.T_START,c.T_END,
                                no_member=c.no_member,
                                )

    era5_var = open_era5(c.DIRECTORY_ERA, c.FILE_ERA, c.VAR)
    era5_var = era5_var.sel(time=slice(sens.time[0].dt.strftime("%Y-%m"),
                                       sens.time[-1].dt.strftime("%Y-%m")))
    era5_var['time']=ctl.time
        
    """
get variables
    """
    ctl = ctl.get(c.VAR)
    sens = sens.get(c.VAR)
    era5_var = era5_var.get(c.VAR)
    
    """
season selection, Clim is the yearly mean
    """
    ctl = dataset_season(ctl, season)
    sens = dataset_season(sens, season)
    era5_var = dataset_season(era5_var, season)
    
    return ctl, sens, era5_var

###############################################################################
###FUNCTIONS FOR DATA MANIPULATION###
###############################################################################
def dataset_season(dset,season):
    """
    SEASON SELECTION
    """
    if season == 'DJF':
        mon = [12,1,2]
    elif season == 'JJA':
        mon = [6,7,8]
    elif season == 'Clim':
        mon = [11,12,1,2,3,4,5,6,7,8,9,10]
        
    dset = dset.sel(time=dset.time.dt.month.isin(mon))
    dset = dset.rolling(time = len(mon)).mean('time')
    dset = dset.sel(time = dset.time.dt.month == mon[-1])
    
    return dset

###############################################################################
def exp_concat(dset1,dset2,dim,lead_exp,season):
    """
    concatenates two datasets
    """

    dset = xr.concat([dset1, dset2], dim=dim)
    dset.to_netcdf(c.RUN_DIR+c.NAME_CTL+'_'+c.VAR+'_lead'+str(lead_exp)+'_'+season+'_s'+str(c.T_START)+'_concat.nc')
    
###############################################################################
def align_time(dset,ref):
    """
    align time
    """
    
    t_ref=ref['time'][:len(dset.time)]
    ref = ref.sel(time=t_ref)
    ref['time']=dset.time
    return ref

###############################################################################
def cross_section(dset):
    """
    ensemble mean
    """
#    if 'member' in dset.dims:
#        dset = dset.mean('member')    
    
    """
    time mean
    """
    #dset = dset.mean('time')
    
    """
    lat average
    """
    LatIndexer = 'lat'
    weights = np.cos(np.deg2rad(dset.lat))
    weights.name = "weights"
    
    dset = dset.sel(**{LatIndexer: slice(-5,5)}).weighted(weights).mean(LatIndexer)
    
    return dset

###############################################################################
def convert_longitudes(ds, lon_name):
# Adjust lon values to make sure they are within (-180, 180)
    ds['_longitude_adjusted'] = xr.where(
        ds[lon_name] > 180,
        ds[lon_name] - 360,
        ds[lon_name])

# reassign the new coords to as the main lon coords
# and sort DataArray using new coordinate values
    ds = (
        ds
        .swap_dims({lon_name: '_longitude_adjusted'})
        .sel(**{'_longitude_adjusted': sorted(ds._longitude_adjusted)})
        .drop(lon_name))

    ds = ds.rename({'_longitude_adjusted': lon_name})
    return ds


###############################################################################
def mse(dataset, ref):
     
    """mse"""
    delta = dataset - ref

    mean_square_error = (delta*delta).mean('time')
    
    return mean_square_error
