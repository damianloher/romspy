#! /usr/bin/env python

import calendar
import netCDF4
import numpy as np
import romspy.Global_settings as Global_settings

ROMSPY_INPUT = Global_settings.Ocean_3d_input
MASK_FILE = ROMSPY_INPUT + '/GCB_RunB_2024_gfort.pop.h.0001.nc'
# File containing the DIC data to be masked:
DIC_FILE = ROMSPY_INPUT + '/CESM_DIC_ann_1751-2023_1x1.nc'

# Read mask for Baltic, Black, Caspian and Red Sea (on POP grid):
nc = netCDF4.Dataset(MASK_FILE,'r')
#spco2 = nc.variables['spco2']
pco2dat = nc.variables['pCO2SURF'][0,:]
mask_2d = (pco2dat.data==0)  # shape yx
nc.close()

# Create 3D mask:
nc = netCDF4.Dataset(DIC_FILE,'r')
nz = len(nc.dimensions['z_t'])
ny = len(nc.dimensions['nlat'])
nx = len(nc.dimensions['nlon'])
nt = len(nc.dimensions['time'])
mask_3d = np.zeros((nz,ny,nx), dtype=bool)
for z in range(nz):
    mask_3d[z,:] = mask_2d
nc.close()

# Process DIC file:
print('Mask DIC in '+DIC_FILE)
nc = netCDF4.Dataset(DIC_FILE,'a')
fill_val = nc.variables['DIC'].missing_value
vobj = nc.variables['DIC']
for t in range(nt):
    dat = nc.variables['DIC'][t,:]
    dat[mask_3d] = fill_val
    vobj[t,:] = dat
# Time variable:
print('Fix time variable')
ndays = 0
tval = []
for yr in range(1751,2024):
    if calendar.isleap(yr):
        tval.append(ndays+366/2)
        ndays += 366
    else:
        tval.append(ndays+365/2)
        ndays += 365
vobj = nc.variables['time']
vobj[:] = tval
vobj.units = 'days since 1751-01-01 00:00:00'
vobj.calendar = 'standard'
nc.close()

print('done')