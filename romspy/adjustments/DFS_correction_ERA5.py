import netCDF4
import cdo
import os
import numpy as np
from romspy import Global_settings

"""
Implements the Drakkar correction based on the Romstools code.
"""

def regrid_dfs_to_romsgrid(grd_file,out_dir,cdo_options,**kwargs):
    """
    Regrids the DFS correction factors to the ROMS grid
      Also sets the correction factors to 1 south of 60S.
    :param grd_file: path to ROMS grid file
    :param out_dir: output folder
    :param cdo_options cdo options
    Can have any of the following optional arguments:
        verbose - whether to print runtime information - default false
    """
    # Check if output file exists: assume nothing needs to be done
    # if this is the case:
    roms_setup = kwargs.get('ROMS_setup', 'pactcs30')
    out_file = f'{out_dir}/{roms_setup}_dfs_factors_roms.nc'
    if os.path.exists(out_file):
        print('found DFS correction factors on ROMS grid: '+out_file)
        return out_file
    verbose = kwargs.get('verbose', False)
    if verbose:
        print('interpolate DFS correction factors to ROMS grid')
    # Create output folder if it does not exist:
    if not os.path.exists(out_dir):
        os.makedirs(out_dir)
    # Shortwave correction factors:
    file_sw = Global_settings.DFS_data_dir + 'SW_dfs_factor.nc'
    file_lw = Global_settings.DFS_data_dir + 'LW_dfs_factor.nc'
    lcdo = cdo.Cdo(debug=verbose)
    outfiles = []
    outfile = lcdo.remapbil(grd_file, input=file_lw, options=cdo_options)
    outfiles.append(outfile)
    outfile = lcdo.remapbil(grd_file, input=file_sw, options=cdo_options)
    outfiles.append(outfile)
    lcdo.merge(input=' '.join(outfiles), output=out_file)
    # Get land mask:
    nc = netCDF4.Dataset(grd_file,'r')
    mask_rho = nc.variables['mask_rho'][:]
    nc.close()
    land_mask = (mask_rho==0)
    # Set correction factors to 1 south of 60N, and set land to 1:
    nc_out = netCDF4.Dataset(out_file,'a')
    lat = nc_out.variables['lat_rho'][:]
    v_dswr = nc_out.variables['factor_dswr']
    v_dlwr = nc_out.variables['factor_dlwr']
    for t in range(len(nc_out.dimensions['time'])):
        factor_dswr = v_dswr[t,:]
        factor_dswr[lat<-60] = 1.0
        factor_dswr[land_mask] = 1.0
        factor_dlwr = v_dlwr[t,:]
        factor_dlwr[lat<-60] = 1.0
        factor_dlwr[land_mask] = 1.0
        v_dswr[t,:] = factor_dswr
        v_dlwr[t,:] = factor_dlwr
    nc_out.close()
    return out_file


def interpolate_clim_dfs_factors_to_daily(output_dir,infile,**kwargs):
    """
    Interpolates the monthly DFS correction factors (assumed to be on the ROMS grid) to daily resolution
    :param output_dir: folder which output file will be read from or written to
    :param grd_file: ROMS grid file
    Can have any of the following optional arguments:
        verbose - whether to print runtime information - default false
        ROMS_setup - name of ROMS setup - default 'pactcs30'
    """
    import scipy.interpolate
    roms_setup = kwargs.get('ROMS_setup', 'pactcs30')
    # Check if the file containing the daily factors is present already:
    out_fname = f'{output_dir}/{roms_setup}_dfs_factors_roms_366_days.nc'
    # Return if this file already exists:
    if os.path.exists(out_fname):
        print('found file with daily DFS correction factors: '+out_fname)
        return out_fname
    verbose = kwargs.get('verbose', False)
    if verbose:
        print('compute daily DFS correction factors from monthly')
    vlist = ['factor_dswr', 'factor_dlwr']
    molen = np.array([31,29,31,30,31,30,31,31,30,31,30,31])
    molen2 = np.zeros((13))
    molen2[1:] = molen
    midmo = np.cumsum(molen2[:-1])+0.5*molen
    # Extend midmo to prevent extrapolation:
    midmo = np.concatenate(([-15.5], midmo, [381.5]))
    mid_daily = np.arange(0.5,366,1)
    # Get dimension sizes:
    nc = netCDF4.Dataset(infile, 'r')
    ny = len(nc.dimensions['eta_rho'])
    nx = len(nc.dimensions['xi_rho'])
    nc.close()
    # Create output file, dimensions and variables:
    if verbose:
        print('creating '+out_fname)
    nc_out = netCDF4.Dataset(out_fname,'w')
    nc_out.createDimension('xi_rho', size=nx)
    nc_out.createDimension('eta_rho', size=ny)
    nc_out.createDimension('time', size=366)
    # time variable:
    vobj = nc_out.createVariable('time', 'f4', ('time',))
    vobj.cycle_length = 366.0
    vobj.calendar = "365_day"
    vobj[:] = mid_daily
    # factor_dlwr:
    vobj = nc_out.createVariable('factor_dlwr', 'f4', ('time','eta_rho', 'xi_rho'), fill_value=-1e20)
    vobj.long_name = "Downward longwave radiation correction factor"
    vobj.description = "Correction factors derived from Drakkar Forcing Set 5.2"
    # factor_dswr:
    vobj = nc_out.createVariable('factor_dswr', 'f4', ('time','eta_rho', 'xi_rho'), fill_value=-1e20)
    vobj.long_name = "Downward shortwave radiation correction factor"
    vobj.description = "Correction factors derived from Drakkar Forcing Set 5.2"
    # Global attributes:
    import datetime
    now = datetime.datetime.today()
    nc_out.remark1 = f'Created using ROMSpy on {now.strftime("%d-%b-%Y %H:%M")}'
    import subprocess
    result = subprocess.run(['git', 'rev-list', '--max-count=1', 'HEAD', 'HEAD'],
                            cwd=os.path.dirname(__file__), stdout=subprocess.PIPE)
    git_commit = result.stdout[:-1].decode('utf-8')
    nc_out.remark2 = 'ROMSpy commit: '+git_commit
    # Do the interpolations:
    for var in vlist:
        if verbose:
            print('   interpolate '+var+'...')
        nc_in = netCDF4.Dataset(infile,'r')
        vobj_in = nc_in.variables[var]
        vobj_out = nc_out.variables[var]
        for j in range(ny):
            for i in range(nx):
                tmp = vobj_in[:,j,i]
                tmp = np.concatenate(([tmp[-1]], tmp, [tmp[0]]))
                func = scipy.interpolate.interp1d(midmo, tmp)
                vobj_out[:,j,i] = func(mid_daily)
        nc_in.close()
    nc_out.close()
    return out_fname


def create_era_frc_radiation(era_path,grd_file,weight_file,out_dir,year,month,
                             era_nr_trecs_day,timavg,cdo_options,**kwargs):
    """
    Creates the ERA radiation data file for a specified year and month
    :param era_path: path to ERA data files
    :param grd_file: path to ROMS grid file
    :param out_dir: output folder
    :param year: year for which to create the data
    :param month: month for which to create the data
    :param era_nr_tsteps_day: number of time records per day in ERA data
    :param timavg [integer]: time resolution of output data. Time averages are computed over blocks of timavg time records
       of the ERA data. E.g. if the ERA5 input files are hourly, the timavg determines
       the time resolution of the output data in hours. If the ERA5 input files are days, timavg determines
       the time resolution of the output data in days.
    :param cdo_options: string of options passed to cdo
    Can have any of the following optional arguments:
        verbose - whether to print runtime information - default false
    """
    verbose = kwargs.get('verbose', False)
    if verbose:
        print(f'create ERA5 radiation files for year {year}, month {month}')
    lcdo = cdo.Cdo(debug=verbose)
    # read grid information:
    nc = netCDF4.Dataset(grd_file,'r')
    Mp = len(nc.dimensions['eta_rho'])
    Lp = len(nc.dimensions['xi_rho'])
    nc.close()
    # Fill value:
    fval = 9.9692099683868690e+36
    # Determine output file name:
    grdf = grd_file.split('/')[-1]
    roms_setup = grdf.split('_')[0]
    out_fname = f'{out_dir}/{roms_setup}_ERA5_radiation_{year}_{month:02}.nc'
    # Create output file, overwriting any existing file with the
    # same name:
    if verbose:
        print(f'create ERA5 frc radiation file: {out_fname}')
    nc_out = netCDF4.Dataset(out_fname,'w')
    nc_out.createDimension('xi_u', size=Lp-1)
    nc_out.createDimension('eta_u', size=Mp)
    nc_out.createDimension('xi_v', size=Lp)
    nc_out.createDimension('eta_v', size=Mp-1)
    nc_out.createDimension('xi_rho', size=Lp)
    nc_out.createDimension('eta_rho', size=Mp)
    nc_out.createDimension('dlw_time', size=None)
    nc_out.createDimension('sthr_time', size=None)
    # Create time variables:
    for var in ['dlw_time','sthr_time']:
        vobj = nc_out.createVariable(var, 'f4', (var,))
        #vobj.units = "days since {}-01-01 00:00".format(yr)
        vobj.climatological = "means"
    # Create data variables:
    vobj = nc_out.createVariable('dlwrad', 'f4', ('dlw_time','eta_rho', 'xi_rho'), fill_value=fval)
    vobj.long_name = 'downward thermal longwave radiation'
    vobj.units = 'Watt meter-2'
    vobj.missing_value = fval
    vobj = nc_out.createVariable('sthrad', 'f4', ('sthr_time','eta_rho', 'xi_rho'), fill_value=fval)
    vobj.long_name = 'net thermal radiation'
    vobj.units = 'Watt meter-2'
    vobj.missing_value = fval
    # Interpolate to the ROMS grid:
    #   str: surface net thermal radiation
    #   strd: surface longwave downward radiation
    for var in ['str','strd']:
        if var == 'str':
            vobj_out = nc_out.variables['sthrad']
            vobj_out_time = nc_out.variables['sthr_time']
            infile = f'{era_path}/{year}/ERA5_{year}_{month:02}.nc'
            if not os.path.exists(infile):
                infile = f'{era_path}/for_romstools/{year}/ERA5_{year}_{month:02}.nc'
                if not os.path.exists(infile):
                    infile = f'{era_path}/{year}/ERA5_{year}_{month:02}_daily.nc'
                    if not os.path.exists(infile):
                        msg = 'no ERA5 data found for variable "str" (surface net thermal radiation)'
                        raise ValueError(msg)
        else:
            vobj_out = nc_out.variables['dlwrad']
            vobj_out_time = nc_out.variables['dlw_time']
            infile = f'{era_path}/{year}/ERA5_{year}_{month:02}_strd.nc'
            if not os.path.exists(infile):
                infile = f'{era_path}/for_romstools/{year}/ERA5_{year}_{month:02}_strd.nc'
                if not os.path.exists(infile):
                    infile = f'{era_path}/{year}/ERA5_{year}_{month:02}_strd_daily.nc'
                    if not os.path.exists(infile):
                        infile = f'{era_path}/{year}/ERA5_{year}_{month:02}.nc'
                        if not os.path.exists(infile):
                            infile = f'{era_path}/for_romstools/{year}/ERA5_{year}_{month:02}.nc'
                            if not os.path.exists(infile):
                                infile = f'{era_path}/{year}/ERA5_{year}_{month:02}_daily.nc'
                                if not os.path.exists(infile):
                                    msg = 'no ERA5 data found for variable "strd" (surface longwave downward radiation)'
                                    raise ValueError(msg)
        outfile = lcdo.remap(grd_file + ',' + weight_file, input=(f' -selname,{var}' + ' ' + infile),
                    options=cdo_options)
        # Scale variable in outfile:
        nc_in = netCDF4.Dataset(outfile, 'r')
        nt = len(nc_in.dimensions['time'])
        #tmp = nc_in.variables[var][:]
        vobj_in = nc_in.variables[var]
        vobj_in_time = nc_in.variables[vobj_in.dimensions[0]]
        # Conversion factor to convert to time averages [sec^-1]:
        stimei = era_nr_trecs_day/(24*3600.0)
        if vobj_in_time.units[:5] == 'hours':
            #vobj_out_time[:] = vobj_in_time[:]/24.0
            tscale = 24.0  # scaling factor to convert time values to days
            vobj_out_time.units = 'days since ' + vobj_in_time.units.split(' ')[-1]
        elif vobj_in_time.units[:4] == 'days':
            #vobj_out_time[:] = vobj_in_time[:]
            tscale = 1.0
            vobj_out_time.units = vobj_in_time.units
        elif vobj_in_time.units.startswith("seconds"):
            tscale = 86400.0
            vobj_out_time.units = 'days since ' + vobj_in_time.units.split(' ')[-1]
        else:
            msg = f'time unit of variable "{vobj_in_time.name}" not supported: {vobj_in_time.units}'
            raise ValueError(msg)
        # Average over time records and scale data:
        t1 = 0
        t2 = timavg
        tidx = 0
        while t2 <= nt:
            # Convert time values to days:
            vobj_out_time[tidx] = np.mean(vobj_in_time[t1:t2])/tscale
            vobj_out[tidx,:] = stimei*np.mean(nc_in.variables[var][t1:t2,:], axis=0)
            t1 = t2
            t2 += timavg
            tidx += 1
        nc_in.close()
        # Remove temporary file produced by cdo:
        os.system(f'rm -f {outfile}')
    # Close output file:
    nc_out.close()
    if verbose:
        print('created radiation file: '+out_fname)
    # Return path to output file:
    return out_fname


def make_drakkar_correction(file,DFS_corr_factor_file,rad_file,month,tsteps_perday,**kwargs):
    """
    Performs the Drakkar correction on a single file
    :param file: forcing file on which the correction is to be done
    :param DFS_corr_factor_file: file containing the correction factor on the ROMS grid
    :param rad_file: file containing ERA5 radiation data
    :param month: month covered by forcing file given by argument "file" (1=January, 12=December)
    :param tsteps_perday: number of time records per day in forcing file given by argument "file"
    Can have any of the following optional arguments:
        verbose - whether to print runtime information - default false
    """
    # Check first if the Drakkar correction can be done at all:
    nc_frc = netCDF4.Dataset(file,'r')
    if not ("swrad" in nc_frc.variables and "shflux" in nc_frc.variables):
        print(f'   skip Drakkar correction: one or both of swrad and shflux not present in {file}')
        nc_frc.close()
        return
    nc_frc.close()
    verbose = kwargs.get('verbose', False)
    if verbose:
        print('apply Drakkar correction to '+file)
    # Determine the first time record to use for the correction
    # factors (assuming they are daily):
    molen = np.array([0,31,29,31,30,31,30,31,31,30,31,30,31])
    nrdays = np.cumsum(molen)
    t_dfs_ini = nrdays[month-1]
    # Open Netcdf files and loop over the forcing time records:
    nc_frc = netCDF4.Dataset(file,'a')
    nc_corr_fac = netCDF4.Dataset(DFS_corr_factor_file,'r')
    nc_rad = netCDF4.Dataset(rad_file,'r')
    nt = len(nc_frc.dimensions['time'])
    for t in range(nt):
        ssr = nc_frc.variables['swrad'][t,:]
        str = nc_rad.variables['sthrad'][t,:]
        strd = nc_rad.variables['dlwrad'][t,:]
        shflux = nc_frc.variables['shflux'][t,:]
        t_dfs = t_dfs_ini + t//tsteps_perday
        #if verbose:
        #    print(f'   t_ROMS_frc = {t}, t_DFS = {t_dfs}')
        factor_swrad = nc_corr_fac['factor_dswr'][t_dfs,:]
        factor_dlwrad = nc_corr_fac['factor_dlwr'][t_dfs,:]
        # Do the correction:
        ssr_corr = ssr*factor_swrad
        # Compute sensible and latent heat fluxes:
        sens_lat = shflux - ssr - str
        # Thermal radiation upward from uncorrected terms:
        stru = strd - str
        # Correct longwave downward flux:
        strd_corr = factor_dlwrad*strd
        # Correct net thermal radiation:
        str_corr = strd_corr - stru
        # Correct net heat flux:
        shflux_corr = sens_lat + ssr_corr + str_corr
        # Write corrected fluxes to forcing file:
        nc_frc['shflux'][t,:] = shflux_corr
        nc_frc['swrad'][t,:] = ssr_corr
    nc_corr_fac.close()
    nc_rad.close()
    nc_frc.close()
