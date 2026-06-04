"""
Author: Nicolas Munnich
License: GNU GPL2+
"""

import netCDF4
from romspy import UP_data_paths
import sys
import os


def str_adjustment(file: str, group_files: str, cdo, options, verbose, **kwargs):
    """
    Add attributes to sustr and svstr and correct unit (if necessary)
    :param preprocessor: PreProcessor object, contains a lot of necessary adjustments
    :param file: file which contains all necessary parameters to calculate dQdSST
    :param group_files: list of all files which contain necessary parameters to calculate dQdSST
    :param flags: any extra info necessary
    :return:
    """
    if verbose:
        print('Adjusting for surface stress')
    # Conversion factor N/m^2 s --> N/m^2:
    scale = kwargs['wind_stress_scale']
    temp = cdo.setattribute(
        "sustr@long_name='surface u-momentum stress',svstr@long_name='surface v-momentum stress',"
        "sustr@units='Newton meter-2',svstr@units='Newton meter-2',sustr@files='"
        + group_files.replace(',',' ') + "',svstr@files='" + group_files.replace(',',' ') + "'",
        input=f"-aexpr,'sustr={scale}*ewss;svstr={scale}*nsss;' {file}", options=options)
    cdo.copy(input=temp, output=file, options=options)

    if verbose:
        print('surface stress has been adjusted')


def dust_adjustment(file: str, group_files: str, cdo, options, verbose, **kwargs):
    """
    Calculate iron
    :param preprocessor: PreProcessor object, contains a lot of necessary adjustments
    :param file: file which contains all necessary parameters to calculate dQdSST
    :param group_files: list of all files which contain necessary parameters to calculate dQdSST
    :param flags: any extra info necessary
    :return:
    """
    if verbose:
        print('Adjusting dust to add iron')

    iron_factor = 62668.0
    temp = cdo.setattribute(
        "iron@long_name='iron_flux',iron@units='nmol/cm2/s',iron@files='" + group_files + "'",
        input="-aexpr,'iron=" + str(iron_factor) + " * dust' " + file,
        options=options)
    cdo.copy(input=temp, output=file, options=options)

    if verbose:
        print('Iron has been added')


def swrad_adjustment(file: str, group_files: str, cdo, options, verbose, **kwargs):
    """
    Scale solar flux (to convert from J/m^2 to W/m^2, which is expected by ROMS)
    :param preprocessor: PreProcessor object, contains a lot of necessary adjustments
    :param file: file which contains all necessary parameters to calculate dQdSST
    :param group_files: list of all files which contain necessary parameters to calculate dQdSST
    :param flags: any extra info necessary
    :return:
    """
    if verbose:
        print("Adjusting for swrad")
    scale = kwargs['swrad_scale']
    temp = cdo.aexpr(f"swrad={scale}*ssr", input="-setattribute,swrad@units='Watt meter-2' "+file, options=options)
    cdo.copy(input=temp, output=file, options=options)
    if verbose:
        print('swrad scaled')

def swflux_adjustment(file: str, group_files: str, cdo, options, verbose, **kwargs):
    """
    Calculate surface fresh water flux
    :param preprocessor: PreProcessor object, contains a lot of necessary adjustments
    :param file: file which contains all necessary parameters to calculate dQdSST
    :param group_files: list of all files which contain necessary parameters to calculate dQdSST
    :param flags: any extra info necessary
    :return:
    """
    if verbose:
        print("Adjusting for swflux")
        sys.stdout.flush()
    # Get conversion factor [m]/"Integration time" --> cm/day (as a string):
    scale = kwargs['swflux_scale']
    temp = cdo.setattribute(
        "swflux@long_name='surface fresh water flux',swflux@units='cm day-1',swflux@files='" + group_files.replace(',',' ') + "'",
        #input="-delname,evap" + (
        #    # 'flags["include_precip"]' instead of 'flags.get("include_precip", False)' makes precip flag mandatory
        #    ",precip" if not include_precip else "") + " -aexpr,'swflux=" + scale + "*(precip+evap)' " + file,
        input=" -aexpr,'swflux=" + scale + "*(tp+e)' " + file,
        options=options)
    cdo.copy(input=temp, output=file, options=options)

    if verbose:
        print("swflux has been added")
        sys.stdout.flush()


def shflux_adjustment(file: str, group_files: str, cdo, options, verbose, **kwargs):
    """
    Compute the surface net heat flux
    :param preprocessor: PreProcessor object, contains a lot of necessary adjustments
    :param file: file which contains all necessary parameters to calculate dQdSST
    :param group_files: list of all files which contain necessary parameters to calculate dQdSST
    :param flags: any extra info necessary
    :return:
    """
    if verbose:
        print('Adjusting for shflux')
        sys.stdout.flush()
    scale = kwargs['swrad_scale']
    # Check for some unnecessary variables in the output file:
    vlist = ['lwrad','sensheat','latheat']
    vlist_present = []
    nc = netCDF4.Dataset(file,'r')
    for v in vlist:
        if v in nc.variables:
            vlist_present.append(v)
    nc.close()
    if len(vlist_present) > 0:
        input_str = '-delname,' + ','.join(vlist_present)
    else:
        input_str = ''
    # aexpr to produce the new variable
    # delname removes unnecessary variables
    # setattribute sets the necessary attributes of new variable
    temp = cdo.setattribute(
        "shflux@long_name='surface net heat flux',shflux@units='Watt meter-2',shflux@files='"
        + group_files.replace(',',' ') + "'",
        input=input_str + f"-aexpr,'shflux=swrad+{scale}*(sshf+slhf+str)' " + file,
        options=options)
    cdo.copy(input=temp, output=file, options=options)

    if verbose:
        print('shflux has been modified')
        sys.stdout.flush()


def dqdsst_adjustment(file: str, group_files: str, cdo, options, verbose, **kwargs):
    """
        Compute the kinematic surface net heat flux sensitivity to the sea surface temperature: dQdSST
        Q_model ~ Q + dQdSST * (T_model - sst)
        dQdSST = - 4 * eps * stef * T^3 - airdens * Cp * CH * U - airdens * CE * L * U * 2353 * ln(10 * q_s / T^2)
        sst: sea surface temperature (Celsius)
        sat: sea level air temperature (Celsius)
        airdens: sea level air density (kg/m^3)
        w3: wind speed (m/s)
        qsea: sea level specific humidity (g/kg)
        dQdSST: kinematic surface net heat flux sensitivity to the sea surface temperature (Watts meter-2 Celsius-1)
        :param preprocessor: PreProcessor object, contains a lot of necessary adjustments
        :param file: file which contains all necessary parameters to calculate dQdSST
        :param group_files: list of all files which contain necessary parameters to calculate dQdSST
        :param flags: any extra info necessary
    """

    if verbose:
        print("Calculating dQdSST")
        sys.stdout.flush()

    # Specific heat of atmosphere
    Cp = str(1004.8)
    # Sensible heat transfer coefficient (stable condition)
    Ch = str(0.66E-3)
    # Latent heat transfer coefficient (stable condition)
    Ce = str(1.15E-3)
    # Emissivity coefficient
    # eps = str(0.98)
    # Stefan constant
    stef = str(5.6697E-8)
    # SST (Kelvin)
    sstk_formula = "_sst = sst + 273.15;"
    # Latent heat of vaporisation (J.Kg-1)
    l = "(2.5008E6 - 2.3E3 * sat)"
    # Infrared contribution
    q1 = "(-4.0 * " + stef + "* _sst * _sst * _sst)"
    # Sensible heat contribution
    q2 = "(-airdens * " + Cp + " * " + Ch + "* w3)"
    # Latent heat contribution
    # 2.30258509299404590109361379290930926799774169921875 = ln(10)
    # 5417.982723814990094979293644428253173828125 = ln(10) * 2353
    # Factor 0.001 converts humidity g/kg to kg/kg:
    dqsdt = "(5417.9827238149900949792936444282531738281250 * 0.001 * qsea / (_sst * _sst))"
    q3 = "(-airdens * " + Ce + " * " + l + " * w3 * " + dqsdt + " )"

    dqdsst_formula = "'" + sstk_formula + "dQdSST=" + q1 + " + " + q2 + " + " + q3 + "'"

    # Check for some unnecessary variables in the output file:
    vlist = ['sat','rho_atm','U','qsea','airdens','w3']
    vlist_present = []
    nc = netCDF4.Dataset(file,'r')
    for v in vlist:
        if v in nc.variables:
            vlist_present.append(v)
    nc.close()
    if len(vlist_present) > 0:
        input_str = '-delname,' + ','.join(vlist_present)
    else:
        input_str = ''
    temp = cdo.setattribute(
        "dQdSST@long_name='surface net heat flux sensitivity to SST',"
        "dQdSST@units='Watts meter-2 Celsius-1',dQdSST@files='" + group_files.replace(',',' ') + "'",
        input=input_str + " -aexpr," + dqdsst_formula + " " + file,
        options=options)
    cdo.copy(input=temp, output=file, options=options)

    if verbose:
        print("dQdSST has been added.")
        sys.stdout.flush()


def coads05_time_axes_adjustment(file: str, group_files: str, cdo, options, verbose, **kwargs):
    """
        Makes sure the time axes and variables in COADS05 based files are as expected by ROMS
        :param preprocessor: PreProcessor object, contains a lot of necessary adjustments
        :param file: file which contains all necessary parameters to calculate dQdSST
        :param group_files: list of all files which contain necessary parameters to calculate dQdSST
        :param flags: any extra info necessary
    """
    sources = kwargs['sources']
    group_index = kwargs['group_index']
    data_source = sources[group_index]['data_source']
    if data_source != 'COADS05':
        print('   skipped: data source is not COADS05')
        return
    import os
    # Check if the current file looks like a COADS05 ROMS forcing file already:
    rename_hdims = False
    taxes = ['sss_time','T']
    nc = netCDF4.Dataset(file,'r')
    # if not 'SSS' in nc.variables:
    #     msg = 'format of forcing file not recoginzed: neither sustr nor SSS found'
    #     raise ValueError(msg)
    taxis_missing = False
    for taxis in taxes:
        if not taxis in nc.dimensions:
            taxis_missing = True
    if 'x' in nc.dimensions or 'y' in nc.dimensions:
        rename_hdims = True
    if not (taxis_missing or rename_hdims):
        # File is assumed to be in a ROMS compatible format, so nothing
        # to do:
        nc.close()
        return
    # Search for a time dimension in the current file:
    tdims = []
    if 'time' in nc.dimensions:
        tdims.append('time')
    elif 'T' in nc.dimensions:
        tdims.append('T')
    if len(tdims) == 0:
        nc.close()
        msg = 'no time dimension found in file: '+file
        raise ValueError(msg)
    elif len(tdims) > 1:
        nc.close()
        msg = f'multiple time dimensions found in file: {file}\n'
        msg += 'Exactly 1 time dimension must be present.'
        raise ValueError(msg)
    nc.close()
    # Rename current file:
    file_old = file + '_old'
    os.system(f'mv {file} {file_old}')
    if verbose:
        print(f'   rename: {file}  -->  {file_old}')
    nc_old = netCDF4.Dataset(file_old,'r')
    nt = len(nc_old.dimensions[tdims[0]])
    if 'y' in nc_old.dimensions:
        ny = len(nc_old.dimensions['y'])
    elif 'eta_rho' in nc_old.dimensions:
        ny = len(nc_old.dimensions['eta_rho'])
    else:
        msg = 'no y dimension found'
        raise ValueError(msg)
    if 'x' in nc_old.dimensions:
        nx = len(nc_old.dimensions['x'])
    elif 'xi_rho' in nc_old.dimensions:
        nx = len(nc_old.dimensions['xi_rho'])
    else:
        msg = 'no x dimension found'
        raise ValueError(msg)
    # Create new file with the same name as the original one:
    nc = netCDF4.Dataset(file,'w')
    # Create dimensions:
    nc.createDimension('xi_rho', size=nx)
    if 'xi_u' in nc_old.dimensions:
        nc.createDimension('xi_u', size=nx-1)
    nc.createDimension('eta_rho', size=ny)
    if 'eta_v' in nc_old.dimensions:
        nc.createDimension('eta_v', size=ny-1)
    for tax in taxes:
        nc.createDimension(tax, size=nt)
    if 'bnds' in nc_old.dimensions:
        nc.createDimension('bnds', size=len(nc_old.dimensions['bnds']))
    if 'nv4' in nc_old.dimensions:
        nc.createDimension('nv4', size=len(nc_old.dimensions['nv4']))
    # Create time variables:
    for tax in taxes:
        vobj = nc.createVariable(tax, 'f4', (tax,))
        vobj.long_name = tax
        vobj.climatological = "means"
        if nt == 12:
            # Assume monthly climatology:
            vobj.units = "days since 1979-01-01 00:00"
            vobj.cycle_length = 365.25
            vobj[:] = [15.5, 45, 74.5, 105, 135.5, 166, 196.5, 227.5, 258, 288.5, 319, 349.5]
        else:
            msg = f'size of dimension "{tax}" must be 12 (found {nt})'
            raise ValueError(msg)
                
    # Create forcing variables with correct time dimension:
    t_map = {'SSS': 'sss_time', 'dQdSST': 'sss_time'}
    for var in nc_old.variables:
        if not var in t_map:
            continue
        # SST should be skipped in the COADS files:
        if var == 'SST':
            continue
        print('   var = '+var)
        vobj_old = nc_old.variables[var]
        fill_val = vobj_old.missing_value
        dims_old = nc_old.variables[var].dimensions
        if tdims[0] in dims_old:
            dims = (t_map[var],)
        for d in dims_old[1:]:
            if d == 'x':
                dims += ('xi_rho',)
            elif d == 'y':
                dims += ('eta_rho',)
            else:
                dims += (d,)
        vobj = nc.createVariable(var, 'f4', dims, fill_value=fill_val, zlib=True,
                                 complevel=1)
        v_attrs = {x: vobj_old.getncattr(x) for x in vobj_old.ncattrs() if x != "_FillValue"}
        vobj.setncatts(v_attrs)
        vobj[:] = vobj_old[:]
    # Global attributes:
    import datetime
    now = datetime.datetime.today()
    nc.remark1 = f'Created using ROMSpy on {now.strftime("%d-%b-%Y %H:%M")}'
    import subprocess
    result = subprocess.run(['git', 'rev-list', '--max-count=1', 'HEAD', 'HEAD'],
                            cwd=os.path.dirname(__file__), stdout=subprocess.PIPE)
    git_commit = result.stdout[:-1].decode('utf-8')
    nc.remark2 = 'ROMSpy commit: '+git_commit
    att_dict = dict()
    # Copy global attributes from old file:
    for att in nc_old.ncattrs():
        att_dict[att] = getattr(nc_old,att)
    nc.setncatts(att_dict)
    nc_old.close()
    nc.close()
    # Remove superfluous files:
    os.system(f'rm -f {file_old} *_humidity_*.nc *_rho_air_*.nc *_t_air_*.nc *_u_air_*.nc')
    os.system('rm -f *_SSS_*.nc *_SST_*.nc')
    if verbose:
        print(f'   wrote {file}')
        sys.stdout.flush()

def era5_extrapolate_adjustment(file: str, group_files: str, cdo, options, verbose, **kwargs):
    sources = kwargs['sources']
    group_index = kwargs['group_index']
    data_source = sources[group_index]['data_source']
    if data_source != 'ERA5':
        print('   skipped: data source is not ERA5')
        return
    grd_file = kwargs.get('target_grid', None)
    if not grd_file:
        print('   skipped: no ROMS grid file given')
        return
    out_dir = sources[group_index].get('auxiliary_folder',None)
    roms_setup = sources[group_index]['ROMS_setup']
    # List of variables to extrapolate:
    vlist = ["sustr","svstr","shflux","swrad","swflux"]
    sys.path.append(os.path.dirname(__file__)+"/../../c++/build")
    import roms_forcing_utils
    roms_frc = roms_forcing_utils.ROMS_frc(roms_setup, grd_file, out_dir, True)
    roms_frc.extrapolate_land(file, grd_file, vlist)

def era5_time_axes_adjustment(file: str, group_files: str, cdo, options, verbose, **kwargs):
    """
        Makes sure the time axes and variables in ERA5 based files are as expected by ROMS,
        and merges all the monthly forcing data of a year into one file.
        :param preprocessor: PreProcessor object, contains a lot of necessary adjustments
        :param file: file which contains all necessary parameters to calculate dQdSST
        :param group_files: list of all files which contain necessary parameters to calculate dQdSST
        :param flags: any extra info necessary
    """
    sources = kwargs['sources']
    group_index = kwargs['group_index']
    data_source = sources[group_index]['data_source']
    if data_source != 'ERA5':
        print('   skipped: data source is not ERA5')
        return
    import os, datetime, calendar
    # Determine month (Jan=1, ..., Dec=12) of this file: assume the filename ends with _xxx.nc,
    # where xxx denotes the month index (January of start_year = 1, February of
    # start_year = 2, ..., January of start_year+1 = 13, ..., December of start_year+1 = 24, etc)
    midx = int(file[-6:-3])
    if midx%12 == 0:
        # It is a December:
        month = 12
    else:
        month = midx%12
    # Return if this is not a December:
    if month != 12:
        print('   skipped: current month is not a December')
        return
    # Dtermine all files containing ERA5 data for this year:
    flist = []
    for m in range(midx-11,midx+1):
        f = file.replace(f'{midx:03}', f'{m:03}')
        flist.append(f)
    # First and last year for which the forcing is to be produced:
    start_year = sources[group_index].get('start_year', None)
    end_year = sources[group_index].get('end_year', None)
    if not start_year or not end_year:
        msg = 'start_year and end_year must be specified'
        raise ValueError(msg)
    # Get year at which the time axes are to start: set it to start_year if
    # the user has not specified it
    start_year_taxes = sources[group_index].get('start_year', None)
    if not start_year_taxes:
        start_year_taxes = start_year
    # Get first and last years ot the whole simulation: set them to start_year
    # and last_year if the user did not specify them
    start_year_run = sources[group_index].get('start_year_run', None)
    if not start_year_run:
        start_year_run = start_year
    end_year_run = sources[group_index].get('end_year_run', None)
    if not end_year_run:
        end_year_run = end_year
    # Determine year of this file:
    year = (midx-1)//12 + start_year
    # Get time resolution in days, if applicable
    tres_str = sources[group_index]['time_resolution']
    if tres_str:
        if tres_str[:2] == '1d':
            tres = 1.0
            trecs_per_day = 1
        elif tres_str[-1] == "h":
            trecs_per_day = 24/int(tres_str[:-1])
            tres = float(tres_str[:-1])/24.0
        else:
            msg = f'format of time resolution not recognized: {tres_str}'
            raise ValueError(msg)
    else:
        msg = 'no specification found for the time resolution of the monthly forcing files'
        raise ValueError(msg)
    use_cyclic_time_axes = sources[group_index].get('use_cyclic_time_axes', True)
    if not use_cyclic_time_axes:
        # This option determines if the time axes of all years are extended by 2
        # records (1 just before the start, 1 just after the end of each year) or
        # if the time axes are only extended by 1 time rec for start_year_run and
        # end_year_run (1 record just before the start of the run, 1 record just
        # after the end of the run)
        extend_taxes_all_years = sources[group_index].get('extend_taxes_all_years', True)
    # Do some checks with this file:
    nc = netCDF4.Dataset(file,'r')
    if not 'sustr' in nc.variables:
        msg = 'format of forcing file not recognized: sustr not found'
        raise ValueError(msg)
    # Search for a time dimension in the current file:
    tdims = []
    if 'time' in nc.dimensions:
        tdims.append('time')
    elif 'T' in nc.dimensions:
        tdims.append('T')
    if len(tdims) == 0:
        nc.close()
        msg = 'no time dimension found in file: '+file
        raise ValueError(msg)
    elif len(tdims) > 1:
        nc.close()
        msg = f'multiple time dimensions found in file: {file}\n'
        msg += 'Exactly 1 time dimension must be present.'
        raise ValueError(msg)
    nc.close()
    nc_old = netCDF4.Dataset(file,'r')
    nt = len(nc_old.dimensions[tdims[0]])
    if 'y' in nc_old.dimensions:
        ny = len(nc_old.dimensions['y'])
    elif 'eta_rho' in nc_old.dimensions:
        ny = len(nc_old.dimensions['eta_rho'])
    else:
        msg = 'no y dimension found'
        raise ValueError(msg)
    if 'x' in nc_old.dimensions:
        nx = len(nc_old.dimensions['x'])
    elif 'xi_rho' in nc_old.dimensions:
        nx = len(nc_old.dimensions['xi_rho'])
    else:
        msg = 'no x dimension found'
        raise ValueError(msg)
    #
    # Create file containing the data of the whole current year:
    #
    fname = f'{file[:-9]}_{year}.nc'
    nc = netCDF4.Dataset(fname,'w')
    # Create dimensions:
    nc.createDimension('xi_rho', size=nx)
    if 'xi_u' in nc_old.dimensions:
        nc.createDimension('xi_u', size=nx-1)
    nc.createDimension('eta_rho', size=ny)
    if 'eta_v' in nc_old.dimensions:
        nc.createDimension('eta_v', size=ny-1)
    if 'bnds' in nc_old.dimensions:
        nc.createDimension('bnds', size=len(nc_old.dimensions['bnds']))
    if 'nv4' in nc_old.dimensions:
        nc.createDimension('nv4', size=len(nc_old.dimensions['nv4']))
    nc_old.close()
    # Time axes:
    taxes = ['sms_time', 'shf_time', 'swf_time', 'srf_time', 'sst_time']
    if calendar.isleap(year):
        nt = 366
        days_per_month = [31,29,31,30,31,30,31,31,30,31,30,31]
    else:
        nt = 365
        days_per_month = [31,28,31,30,31,30,31,31,30,31,30,31]
    if not use_cyclic_time_axes:
        if extend_taxes_all_years:
            # Increase size of time axes by 2:
            nt += 2
        else:
            if year == start_year_run or year == end_year_run:
                nt += 1
    for tax in taxes:
        nc.createDimension(tax, size=nt*trecs_per_day)
    # Create time variables:
    for tax in taxes:
        vobj = nc.createVariable(tax, 'f4', (tax,))
        vobj.long_name = tax
        vobj.climatological = "means"
        if use_cyclic_time_axes:
            vobj.units = f"days since {year}-01-01 00:00"
            if calendar.isleap(year):
                vobj.cycle_length = 366.0
            else:
                vobj.cycle_length = 365.0
        else:
            vobj.units = f"days since {start_year_taxes}-01-01 00:00"
            vobj.cycle_length = 0.0
        # Write time variable for each month to output file:
        for m in range(12):
            f = flist[m]
            nc_old = netCDF4.Dataset(f,'r')
            # Check unit of old time variable (days or hours since some date):
            vobj_old = nc_old.variables[tdims[0]]
            unit_old = vobj_old.units
            unit_old_list = unit_old.split(' ')
            start_year_old = int(unit_old_list[2][:4])
            start_month_old = int(unit_old_list[2][5:7])
            start_day_old = int(unit_old_list[2][8:10])
            start_hour_old = int(unit_old_list[3][:2])
            if start_hour_old > 0:
                msg = f'time axes must start at hour 0 of a day (found hour {start_hour_old})'
                raise ValueError(msg)
            # Adapt time axes:
            if use_cyclic_time_axes:
                d1 = datetime.date(year,1,1)
            else:
                d1 = datetime.date(start_year_taxes,1,1)
            d2 = datetime.date(start_year_old,start_month_old,start_day_old)
            delta = (d1-d2).days
            if unit_old_list[0] == 'days':
                tmp = vobj_old[:] - delta
            elif unit_old_list[0] == 'hours':
                tmp = vobj_old[:]/24.0 - delta
            else:
                msg = f'time unit must be "days" or "hours" (found "{unit_old_list[0]}")'
                raise ValueError(msg)
            # Write time values:
            nt_m = tmp.shape[0]  # number of time records in monthly file
            if m == 0:
                # January:
                if use_cyclic_time_axes:
                    t1 = 0
                elif extend_taxes_all_years or year==start_year_run:
                    t1 = 1
                else:
                    t1 = 0
                t2 = t1 + nt_m-1  # last time record to be written
            else:
                t1 = t2+1
                t2 = t1 + nt_m - 1
            vobj[t1:t2+1] = tmp
            if not use_cyclic_time_axes:
                if m == 0 and (extend_taxes_all_years or year==start_year_run):
                    # January:
                    vobj[0] = vobj[1] - tres
                elif m == 11 and (extend_taxes_all_years or year==end_year_run):
                    # December:
                    vobj[-1] = vobj[-2] + tres
            nc_old.close()
                
    # Create forcing variables with correct time dimension:
    t_map = {'sustr': 'sms_time', 'svstr': 'sms_time', 'shflux': 'shf_time', 'swflux': 'swf_time',
             'swrad': 'srf_time', 'SST': 'sst_time'}
    var_map = {'sst': 'SST'}
    for var in nc_old.variables:
        if not var in t_map.keys():
            continue
        print('   var = '+var)
        sys.stdout.flush()
        nc_old = netCDF4.Dataset(flist[0],'r')
        vobj_old = nc_old.variables[var]
        fill_val = vobj_old.missing_value
        dims_old = nc_old.variables[var].dimensions
        if tdims[0] in dims_old:
            dims = (t_map[var],)
        for d in dims_old[1:]:
            if d == 'x':
                dims += ('xi_rho',)
            elif d == 'y':
                dims += ('eta_rho',)
            else:
                dims += (d,)
        vobj = nc.createVariable(var, 'f4', dims, fill_value=fill_val, zlib=True,
                                 complevel=1)
        v_attrs = {x: vobj_old.getncattr(x) for x in vobj_old.ncattrs() if x != "_FillValue"}
        if not use_cyclic_time_axes and "cycle_length" in v_attrs:
            del v_attrs["cycle_length"]
        vobj.setncatts(v_attrs)
        nc_old.close()
        # Write data for each month to output file:
        for m in range(12):
            f = flist[m]
            nc_old = netCDF4.Dataset(f,'r')
            vobj_old = nc_old.variables[var]
            tmp = vobj_old[:]
            # Write values:
            nt_m = tmp.shape[0]  # number of time records in monthly file
            if m == 0:
                # January:
                if use_cyclic_time_axes:
                    t1 = 0
                elif extend_taxes_all_years or year==start_year_run:
                    t1 = 1
                else:
                    t1 = 0
                t2 = t1 + nt_m-1  # last time record to be written
            else:
                t1 = t2+1
                t2 = t1 + nt_m - 1
            vobj[t1:t2+1,:] = tmp
            if not use_cyclic_time_axes:
                if m == 0:
                    # January:
                    if year == start_year_run:
                        # Duplicate the first time record:
                        vobj[0,:] = vobj[1,:]
                    elif extend_taxes_all_years:
                        # Set first time record to the last record from December of the
                        # previous year.
                        # Set last time record of December of previous year to the first 
                        # of this year:
                        tmp1 = int(fname[-7:-3]) - 1
                        tmp1 = str(tmp1)
                        tmp2 = list(fname)
                        tmp2[-7:-3] = tmp1
                        file_prev = "".join(tmp2)
                        if verbose:
                            print('     file of prev Dec: '+file_prev)
                        nc_prev = netCDF4.Dataset(file_prev,'a')
                        vobj_prev = nc_prev.variables[var]
                        vobj_prev[-1,:] = vobj[1,:]
                        vobj[0,:] = vobj_prev[-2,:]
                        nc_prev.close()
                elif m == 11:
                    # December:
                    # We only need to assign the last time record if this is the final
                    # year of the whole run (if necessary, the other years will be
                    # completed when January of the following year is processed)
                    if year == end_year_run:
                        # Duplicate last time record:
                        vobj[-1,:] = vobj[-2,:]
                    elif extend_taxes_all_years and year == end_year:
                        # Try to fill last time record from data of the next year,
                        # or duplicate it, if this data is not found:
                        tmp1 = int(fname[-7:-3]) + 1
                        tmp1 = str(tmp1)
                        tmp2 = list(fname)
                        tmp2[-7:-3] = tmp1
                        file_next = "".join(tmp2)
                        if os.path.exists(file_next):
                            if verbose:
                                print('     file of next Jan: '+file_next)
                            nc_next = netCDF4.Dataset(file_next,'a')
                            vobj_next = nc_next.variables[var]
                            vobj[-1,:] = vobj_next[1,:]
                            nc_next.close()
                        else:
                            # Duplicate last time record:
                            vobj[-1,:] = vobj[-2,:]
            nc_old.close()

    # Global attributes:
    import datetime
    now = datetime.datetime.today()
    nc.remark1 = 'Created using ROMSpy on {}'.format(now.strftime("%d-%b-%Y %H:%M"))
    import subprocess
    result = subprocess.run(['git', 'rev-list', '--max-count=1', 'HEAD', 'HEAD'],
                            cwd=os.path.dirname(__file__), stdout=subprocess.PIPE)
    git_commit = result.stdout[:-1].decode('utf-8')
    nc.remark2 = 'ROMSpy commit: '+git_commit
    # Copy global attributes from January file:
    nc_old = netCDF4.Dataset(flist[0],'r')
    att_dict = dict()
    for att in nc_old.ncattrs():
        att_dict[att] = getattr(nc_old,att)
    nc.setncatts(att_dict)
    nc_old.close()
    nc.close()
    if verbose:
        print(f'wrote {fname}')
        sys.stdout.flush()


def check_sst_unit(file: str, group_files: str, cdo, options, verbose, **kwargs):
    """
        Makes sure SST is in Celsius and not K
        :param preprocessor: PreProcessor object, contains a lot of necessary adjustments
        :param file: file which contains all necessary parameters to calculate dQdSST
        :param group_files: list of all files which contain necessary parameters to calculate dQdSST
        :param flags: any extra info necessary
    """
    nc = netCDF4.Dataset(file,'a')
    if 'SST' in nc.variables:
        vobj = nc.variables['SST']
        if vobj.units == 'K' or vobj.units == 'Kelvin':
            vobj.units = 'Celsius'
            tmp = vobj[:]
            vobj[:] = tmp - 273.15
            if verbose:
                print('corrected SST unit (K --> Celsius) in {}'.format(file))
    nc.close()

def fill_missing(file: str, group_files: str, cdo, options, verbose, **kwargs):
    """
        Extrapolate to land for specififed variables
        :param preprocessor: PreProcessor object, contains a lot of necessary adjustments
        :param file: file which we want to process
        :param group_files: list of all data files
        :param flags: any extra info necessary
    """
    fillvars = kwargs.get('fill_missing', [])
    if isinstance(fillvars,bool) and fillvars:
        fillvars = ['SST','sst','swrad','shflux','sustr','svstr','swflux','SSS','dQdSST']
    # Reduce fillvars to existing variables:
    nc = netCDF4.Dataset(file,'r')
    tmp = fillvars.copy()
    for v in tmp:
        if not v in nc.variables:
            fillvars.remove(v)
    nc.close()
    if len(fillvars) > 0:
        if verbose:
            print('fill_missing: fill {}'.format(fillvars))
        vlist = ','.join(fillvars)
        vlist_2 = [v+'_2' for v in fillvars]
        input1 = '-selvar,' + vlist
        temp = cdo.fillmiss2(input=input1+' '+file, options=options)
        outfile = cdo.merge(input="{} {}".format(temp,file), options=options)
        # The merge command generated x_2 variables for x in fillvars, so remove them:
        cdo.delname(','.join(vlist_2), input=outfile, output=file, options=options)

def Drakkar_correction(file: str, group_files: str, cdo, options, verbose, **kwargs):
    """
    Apply Drakkar correction to a single ERA5 based forcing file
    :param preprocessor: PreProcessor object, contains a lot of necessary adjustments
    :param file: file we want to apply the Drakkar correction to
    :param group_files: list of all files which contain necessary parameters to calculate dQdSST
    :param kwargs: any extra info necessary
    """
    # First check if "file" contains all the variables affected by the Drakkar correction:
    # if this is not the case, then we don't do it
    nc_frc = netCDF4.Dataset(file,'r')
    if not ("swrad" in nc_frc.variables and "shflux" in nc_frc.variables):
        print(f'   skipped: one or both of swrad and shflux not present in {file}')
        nc_frc.close()
        return
    nc_frc.close()
    roms_setup = kwargs['ROMS_setup']
    sources = kwargs['sources']
    group_index = kwargs['group_index']
    data_source = sources[group_index]['data_source']
    if data_source != 'ERA5':
        print('   skipped: data source is not ERA5')
        return
    import os, sys
    # Folder where the Drakkar correction files are stored:
    drakkar_corr_folder = sources[group_index].get('auxiliary_folder',None)
    if not drakkar_corr_folder:
        drakkar_corr_folder = os.path.dirname(file)
    # Determine month (Jan=1, ..., Dec=12) of this file: assume the filename ends with _xxx.nc,
    # where xxx denotes the month index (January of start_year = 1, February of
    # start_year = 2, ..., January of start_year+1 = 13, ..., December of start_year+1 = 24, etc)
    midx = int(file[-6:-3])
    if midx%12 == 0:
        # It is a December:
        month = 12
    else:
        month = midx%12
    # Determine year of this file:
    start_year = sources[group_index].get('start_year', None)
    year = (midx-1)//12 + start_year
    # Check the time resolution of the forcing:
    tres = sources[group_index].get('time_resolution')
    sys.path.append(os.path.dirname(__file__))
    import DFS_correction_ERA5
    if tres == '1d':
        era_path = '/net/sea/work/datasets/gridded/atmosphere/2d/reanalysis/era5/daily'
        # Number of time records to average in input files:
        timavg = 1
        # Number of ROMS forcing time records per day:
        ROMS_frc_tstep_day = 1
    elif tres == '1d_1h':
        era_path = '/net/sea/work/datasets/gridded/atmosphere/2d/reanalysis/era5/hourly'
        timavg = 24
        ROMS_frc_tstep_day = 1
    elif tres[-1] == 'h':
        era_path = '/net/sea/work/datasets/gridded/atmosphere/2d/reanalysis/era5/hourly'
        timavg = int(tres[:-1])
        ROMS_frc_tstep_day = int(24/timavg)
    else:
        raise ValueError(f'time resolution not supported: {tres}')
    # Path to ROMS grid file:
    target_grid = kwargs.get('target_grid', None)
    cdo_options = options
    basedir = os.path.dirname(file)
    weight_file = f'{basedir}/weights/bil_weight_g0.nc'
    outfile = DFS_correction_ERA5.regrid_dfs_to_romsgrid(target_grid,drakkar_corr_folder,
                                                         cdo_options,ROMS_setup=roms_setup,verbose=verbose)
    corr_file = DFS_correction_ERA5.interpolate_clim_dfs_factors_to_daily(drakkar_corr_folder,outfile,
                                                                          ROMS_setup=roms_setup,verbose=verbose)
    rad_file = DFS_correction_ERA5.create_era_frc_radiation(era_path,target_grid,weight_file,basedir,
                                year,month,timavg,cdo_options,verbose=verbose)
    DFS_correction_ERA5.make_drakkar_correction(file,corr_file,rad_file,month,ROMS_frc_tstep_day,
                                                verbose=verbose)
    os.system('rm -f '+rad_file)

def river_swflux_correction(file: str, group_files: str, cdo, options, verbose, **kwargs):
    """
    Correct freshwater flux due to river input
    """
    sources = kwargs['sources']
    group_index = kwargs['group_index']
    data_source = sources[group_index]['data_source']
    if data_source != 'ERA5':
        print('   skipped: data source is not ERA5')
        return
    grd_file = kwargs.get('target_grid', None)
    # Time resolution of ROMS forcing: time_res is in hours
    tres = sources[group_index].get('time_resolution')
    if tres.startswith('1d'):
        time_res = 24
    elif tres.endswith('h'):
        time_res = int(tres[:-1])
    else:
        msg = f"invalid forcing time resolution: {tres}"
        raise ValueError(msg)
    print(f"river_swflux_correction: time_res = {time_res}h")
    # Use the folder for auxiliary files to store the daily river inputs:
    out_dir = sources[group_index].get('auxiliary_folder',None)
    # Determine month of this file:
    midx = int(file[-6:-3])
    if midx%12 != 0:
        # It is not a December:
        print('   skipped: current month is not a December')
        return
    # Check if swflux is present in "file":
    nc = netCDF4.Dataset(file,'r')
    if not "swflux" in nc.variables:
        print(f'   skipped: swflux not present in file {file}')
        nc.close()
        return
    nc.close()
    #import romspy.adjustments.river_runoff
    import os, sys
    sys.path.append(os.path.dirname(__file__)+"/../../c++/build")
    import roms_forcing_utils
    import calendar
    import numpy as np
    import time
    start_year = sources[group_index].get('start_year', None)
    year = (midx-1)//12 + start_year
    # Read in river data and distribute input of each river over a certain subregion
    # of the ROMS domain:
    roms_setup = sources[group_index]['ROMS_setup']
    river_input_file = UP_data_paths.Dai_river_runoff_dir + 'Dai2009_river_discharge_monthly_clim.nc'
    roms_frc = roms_forcing_utils.ROMS_frc(roms_setup, grd_file, out_dir, True)
    t1 = time.time()
    river_input_romsgrid = roms_frc.make_river_freshwater(river_input_file, out_dir, time_res, "gaussian")
    t2 = time.time()
    print(f"time for ROMS_forcing::make_river_freshwater: {t2-t1:.2f} s")
    sys.stdout.flush()
    #outfile = romspy.adjustments.river_runoff.make_river_freshwater(grd_file,roms_setup,out_dir, verbose=verbose)
    # Determine all forcing files of the current year:
    flist = []
    for m in range(midx-11,midx+1):
        f = file.replace(f'_{midx:03}.nc', f'_{m:03}.nc')
        flist.append(f)
    # Read river freshwater input for a leap year:
    print(f"add river swflux to ROMS forcing swflux, using {river_input_romsgrid}")
    nc = netCDF4.Dataset(river_input_romsgrid,'r')
    river_swflux = nc.variables['swflux']
    #if calendar.isleap(year):
    #    swflux = np.concatenate((swflux[:59,:],[0.5*(swflux[58,:]+swflux[59,:])], swflux[59:,:]))
    #nc.close()
    # Add river freshwater input to forcing files of this year:
    t1 = 0
    t2 = 0
    fidx = 1
    # Number of time records per day:
    trecs_per_day = int(24.0 / time_res)
    for f in flist:
        print(f"   File = {f}")
        nc = netCDF4.Dataset(f,'a')
        vobj = nc.variables['swflux']
        if fidx==3 and not calendar.isleap(year):
            # Skip February 29:
            t1 = t2 + trecs_per_day
        else:
            t1 = t2
        t2 = t1 + len(nc.dimensions[vobj.dimensions[0]])
        if verbose:
            print(f'   river swflux[{t1}:{t2},:] used')
            sys.stdout.flush()
        tmp = vobj[:]
        vobj[:] = tmp + river_swflux[t1:t2,:]
        nc.close()
        fidx += 1
    if verbose:
        print('river runoff correction applied to swflux')
        sys.stdout.flush()

def seaice_correction(file: str, group_files: str, cdo, options, verbose, **kwargs):
    sources = kwargs['sources']
    group_index = kwargs['group_index']
    data_source = sources[group_index]['data_source']
    if data_source != 'ERA5':
        print('   skipped: data source is not ERA5')
        return
    aux_dir = sources[group_index].get('auxiliary_folder',None)
    # Determine year of this file:
    midx = int(file[-6:-3])
    if midx%12 != 0:
        # It is not a December:
        print('   skipped: current month is not a December')
        return
    # Time resolution of ROMS forcing: time_res is in hours
    tres = sources[group_index].get('time_resolution')
    if tres.startswith('1d'):
        time_res = 24
    elif tres.endswith('h'):
        time_res = int(tres[:-1])
    else:
        msg = f"invalid forcing time resolution: {tres}"
        raise ValueError(msg)
    print(f"seaice_correction: time_res = {time_res}h")
    sys.stdout.flush()
    # Number of forcing time records per day:
    trecs_per_day = int(24.0 / time_res)
    day_per_trec = time_res / 24.0
    import os
    import numpy as np
    start_year = sources[group_index].get('start_year', None)
    year = (midx-1)//12 + start_year
    roms_setup = sources[group_index]['ROMS_setup'].lower()
    # Get path to daily sea ice data file:
    seaicefile_default = f'{aux_dir}/{roms_setup}_seaice_frc.nc'
    seaicefile = kwargs.get('seaicefile', seaicefile_default)
    if os.path.exists(seaicefile):
        print(f"seaice_correction: sea ice file = {seaicefile}")
        nc_seaice = netCDF4.Dataset(seaicefile,'r')
    else:
        msg = f"seaice correction NOT done: could not find sea ice file: {seaicefile}"
        raise ValueError(msg)
    # Get path to daily snow ice data file:
    snowicefile_default = f'{aux_dir}/{roms_setup}_snowice_frc.nc'
    snowicefile = kwargs.get('snowicefile', snowicefile_default)
    if os.path.exists(snowicefile):
        print(f"seaice_correction: snow ice file = {snowicefile}")
        nc_snowice = netCDF4.Dataset(snowicefile,'r')
    else:
        msg = f"seaice correction NOT done: could not find snow ice file: {snowicefile}"
        raise ValueError(msg)
    # Determine all forcing files of the current year:
    flist = []
    for m in range(midx-11,midx+1):
        f = file.replace(f'{midx:03}', f'{m:03}')
        flist.append(f)
    import roms_forcing_utils
    import time
    grd_file = kwargs.get('target_grid', None)
    out_dir = sources[group_index].get('auxiliary_folder',None)
    roms_frc = roms_forcing_utils.ROMS_frc(roms_setup, grd_file, out_dir, True)
    t1 = time.time()
    roms_frc.make_seaice_correction(flist, year, seaicefile, snowicefile, time_res, True)
    t2 = time.time()
    print(f"time for ROMS_forcing::make_seaice_correction: {t2-t1:.2f} s")


forcing_adjustments = [
    {
        'out_var_names': set(), 'in_var_names': {'sustr', 'svstr'},
        'func': str_adjustment
    },
    {
        'out_var_names': set(), 'in_var_names': {'swrad'}, 'func': swrad_adjustment
    },
    {
        'out_var_names': set(), 'in_var_names': {'swrad', 'shflux'},
        'func': shflux_adjustment
    },
    {
        'out_var_names': {'swflux'}, 'in_var_names': {'evap', 'precip'}, 'func': swflux_adjustment
    },
    {
        'out_var_names': {'iron'}, 'in_var_names': {'dust'}, 'func': dust_adjustment
    },
    {
        'out_var_names': {'dQdSST'}, 'in_var_names': {'sst', 'sat', 'airdens', 'w3', 'qsea'},
        'func': dqdsst_adjustment
    },
    {
        'out_var_names': set(), 'in_var_names': {'SST'}, 'func': check_sst_unit
    },
    {
        'out_var_names': set(), 'in_var_names': set(), 'func': Drakkar_correction
    },
    #{
    #    'out_var_names': set(), 'in_var_names': set(), 'func': fill_missing
    #},
    {
        'out_var_names': set(), 'in_var_names': set(), 'func': river_swflux_correction
    },
    {
        'out_var_names': {'seaice'}, 'in_var_names': {'sustr','svstr','shflux','swflux','swrad'}, 'func': seaice_correction
    },
    {
        'out_var_names': set(), 'in_var_names': set(), 'func': era5_extrapolate_adjustment
    },
    {
        'out_var_names': set(), 'in_var_names': set(), 'func': era5_time_axes_adjustment
    },
    {
        'out_var_names': set(), 'in_var_names': set(), 'func': coads05_time_axes_adjustment
    },
]

#
# dict={
#     'out_var_name':['a'],
#     'in_var_name':['b','c'],
#     'func': b_c_to_a
# }
#
# forcing_adjustments.append(dict)
