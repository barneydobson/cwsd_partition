# -*- coding: utf-8 -*-
"""
Created on Fri Dec  6 15:17:07 2019

@author: bdobson
"""

"""Constants
(note, most of these aren't used.. don't worry!)
"""
M3_S_TO_ML_D = 86.4 
MM_KM2_TO_ML = 1e-3 * 1e6 * 1e3 * 1e-6 # mm->m, km2->m2, m3->l, l->Ml
MM_M2_TO_ML = 1e-3 * 1e3 * 1e-6 # mm->m, m3->l, l->Ml
MM_M2_TO_SIM_VOLUME = MM_M2_TO_ML #SIM volume is by default ML, but can be changed by changing MM_TO_SIM_VOLUME
ML_TO_M3 = 1000
PCT_TO_PROP = 1/100
L_TO_ML = 1e-6
FLOAT_ACCURACY = 1e-7
UNBOUNDED_CAPACITY = 1e10
MAX_INFLOW = 1e10
DT_DAYS = 1
POLLUTANTS = ['misc_pol'] # All assume mg/l
GARDEN_MULTIPLIER = 2 # The amount we assume the 'actual' water demand from gardens to be that must be supplied by both rainfall and mains water
PCT_GARDENS = 0.1 # Percentage of area that is people's gardens
PI = 3.141592653589793
PER_DAY_TO_PER_SECOND = 1/(60*60*24)
PER_DAY_TO_PER_HOUR = 1/24
M3_S_TO_M3_DT = 1 # Assume 1 but can change