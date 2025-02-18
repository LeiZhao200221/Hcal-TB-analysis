import numpy as np
import pandas as pd

# Calibrate TOA from raw to nanoseconds
def toa_calib(column_data,calibration_data,end):
    result = []
    non_zero=0
    for k,this_val in enumerate(column_data):
        if(this_val!=0):
            non_zero+=1
            this_bx = 25*k
            if this_val < calibration_data[0][4]:
                this_bx = 25 * (k + 1)
            result.append(((this_val/1023)*25)+this_bx)
    if(len(result)!=1):
        return 0
    else:
        if(end==1):
            return sum(result)-calibration_data[0][5]
        else:
            return sum(result)

# (Almost) pass-through function for TOT but removes weird events with multiple TOT's or exactly 2048 (which is not yet understood)
def tot_calib(column_data):
    result = []
    non_zero=0
    for this_val in column_data:
        if(this_val!=0):
            non_zero+=1
            if(this_val==2048):
                this_val=0
            result.append(this_val)
    if(len(result)!=1):
        return 0
    else:
        return sum(result)