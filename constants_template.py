#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Nov 29 20:43:24 2022

@author: ema
"""
###constants definitiontion
###exoperiments data
NAME_SENS = 'a52o'
NAME_CTL = 'a1ua'
T_START = 1993
SENS_NUMBER = 10
CTL_NUMBER = 10
VAR =
UNITS =

DIM = "3D"
'surface for 2D variables'
PLEV =

SIGN = 0.95

###calculation information
no_member = [0] #member to exclude

###flag for reference
REF = 

if REF == 'era':
###era5###
    ERA_VAR = 
    DIRECTORY_ERA = 
    FILE_ERA = 

elif REF == 'hadcrut':
###Observations###
    DIRECTORY_OBS = 
    FILE_OBS = 

###Directories
FILE_DIR =

RUN_DIR =

OUT_DIR = 

###############################################################################
###flags for diagnostics###
###############################################################################
#mean bias
M_BIAS = True

###############################################################################
#correlation
CORRELATION = False

##############################################################################
#crossection
XSECT = False

###############################################################################
#mean squerad skill score
MSSS = False

###############################################################################
#bootstrap method
BOOTSTRAP = 'standard'

###############################################################################
###climate indices
###############################################################################
###NAO
NAO = False

### AO
AO = False
