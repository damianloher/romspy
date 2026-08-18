#! /usr/bin/env python

# Create 1h wind stress forcing files for Humpac15 setup.

import os
import sys

#sys.path.append('/home/loher/.local/lib/python3.9/site-packages')
sys.path.append('/home/loher/python/romspy')
sys.path.append('/workspace/code/romspy')
sys.path.append('/workspace/code/romspy/romspy')

from romspy import PreProcessorFrc, forcing_adjustments
from romspy import PreProcessorClm, clim_adjustments
from romspy import PreProcessorIni, ini_adjustments
from romspy import PreProcessorBry, bry_adjustments
from romspy import clim_adjustments
from romspy import Global_settings

# Name of ROMS setup:
ROMSsetup = 'humpac15'
if Global_settings.use_container:
    # Path to ROMS grid file:
    ROMS_grid = f"{Global_settings.ROMS_grd}/humpac15_grd.nc"
    # Folder where output files will be stored:
    outdir = '/workspace/ROMSpy_output'
    # Path to sea and snow ice files:
    seaice_file = f"{Global_settings.ROMS_files}/humpac15_seaice_frc.nc"
    snowice_file = f"{Global_settings.ROMS_files}/humpac15_snowice_frc.nc"
    # ERA5 land-sea mask file:
    era_lsm_file = f"{Global_settings.era_path}/land_seamask.nc"
    # Folder for aux files:
    aux_folder = '/workspace/input/aux_files'
else:
    # Path to ROMS grid file:
    ROMS_grid = "/net/sea/work/jahaerri/roms/inputs/humpac15_Ncycle/grd/humpac15_grd.nc"
    # Folder where output files will be stored:
    outdir = '/net/sea/work/loher/ROMSpy_output/ERA5_hourly_6'
    # Path to sea and snow ice files:
    seaice_file = "/net/sea/work/koehne/roms/inputs/humpac15/hindcast_1979_2019/frc/corr_fields_seaice/365days/humpac15_seaice_frc.nc"
    snowice_file = "/net/sea/work/koehne/roms/inputs/humpac15/hindcast_1979_2019/frc/corr_fields_seaice/365days/humpac15_snowice_frc.nc"
    # ERA5 land-sea mask file:
    era_lsm_file = "/net/sea/work/datasets/gridded/atmosphere/2d/reanalysis/era5/land_seamask.nc"
    # Folder for aux files:
    aux_folder = '/net/sea/work/loher/ROMSpy_output/auxiliary_files'
era_lsm_var = "var172"
era_lsm_limit = 0.1

# The following list of dictionaries contains the information about what data is to
# be used for which variable.
sources = [
    # 1) Setting related to atm forcing:
    # ERA5 variables:
    {
        'ROMS_setup': ROMSsetup,
        'data_source': 'ERA5',
        'variables': [
            {'out': 'sustr', 'in': 'ewss'},
            {'out': 'svstr', 'in': 'nsss'},
            {'out': 'swrad', 'in': 'ssr'},
            {'out': 'shflux', 'in': ["sshf","slhf","str"], 'expr': 'shflux=sshf+slhf+str;'},
            {'out': 'evap', 'in': 'e'},
            {'out': 'precip', 'in': 'tp'},
            {'out': 'SST', 'in': 'sst'}
        ],
        'base_folder': f'{Global_settings.era_path}/',
        'interpolation_method': 'bil',
        'time_resolution': '1h',   # '1d': daily forcing, '4h': 4h resolution
        'start_year_run': 1979,  # defines the start of all forcing time axes
        'end_year_run': 2025,
        'start_year': 1979,  # must be equal to start_year_run!
        'end_year': 1984,
        'use_cyclic_time_axes': False, # whether the time axes in the forcing are cyclic or not
        'extend_taxes_all_years': True,
        'auxiliary_folder': aux_folder,
    },
]

print("Output folder: "+outdir)
os.makedirs(outdir, exist_ok=True)
# Forcing files:
preproc_frc = PreProcessorFrc(ROMSsetup, outdir, ROMS_grid, sources, layers=64, hc=250.0,
                           tcline=250.0, theta_s=10.0, theta_b=4,
                           seaice_file=seaice_file, snowice_file=snowice_file,
                           fillmiss_after_hor=False, use_ROMS_grdfile=True,
                           lsm_file=era_lsm_file, lsm_var=era_lsm_var, lsm_limit=era_lsm_limit,
                           verbose=1)
preproc_frc.adjustments = forcing_adjustments
preproc_frc.mark_as_vectors('sustr','svstr')
preproc_frc.obc = [1, 0, 0, 0] # open at [S, E, N, W]
preproc_frc.make()
