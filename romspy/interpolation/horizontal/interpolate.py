from itertools import count
import netCDF4
import numpy as np
import os
import romspy.Global_settings as Global_settings
import subprocess
import uuid
import time
import sys

"""
Author: Nicolas Munnich
License: GNU GPL2+
"""


def cdo_interpolate(cdo, file: str, weight: str, target_grid: str, variables: list, all_files: str, options: str,
                outdir: str, outfile: str = None, lsm_file: str = "", lsm_var: str = "", lsm_limit: float = 0.0,
                verbose: int = 1) -> str:
    """
    Extract a variable from file and interpolate according to weight
    :param cdo: cdo object
    :param file: file to extract stuff from
    :param weight: precalculated interpolation weight
    :param target_grid: Grid to interpolate onto
    :param variables: list of variables to extract from the file, paired if the name in the file is different
    :param all_files: list of all files used as a source for these variables
    :param options: cdo options
    :param outfile: name of file to output. if none a temporary file is generated.
    :param verbose: whether to print information
    """
    if len(variables) == 0:
        return
    if verbose > 0:
        print('Horizontally interpolating ' + " ".join([x['out'] for x in variables]) + " from " + file)

    import netCDF4
    try:
        nc = netCDF4.Dataset(file,'r')
        nc.close()
        return __interpolate_nc(cdo, file, weight, target_grid, variables, all_files, outfile, outdir,
                        options, lsm_file, lsm_var, lsm_limit, verbose)
    except OSError:
        #extension = "." + file[::-1].split(".")[0][::-1]  # get file extension
        #raise ValueError(extension + " source files are currently not supported. "
        raise ValueError("not a Netcdf file: " + file + "\n" +
                         "If you wish to add this functionality, " +
                         "please look under romspy/interpolation/horizontal/interpolate.py")
    # Room for additional/alternative conversion functions if they appear necessary


def __interpolate_nc(cdo, file: str, weight: list, target: str, variables: list, all_files: str, outfile: str,
                     outdir: str, options: str, lsm_file: str, lsm_var: str, lsm_limit: float, verbose: int) -> str:
    """
    Interpolate and extract variable
    :param cdo: cdo object
    :param file: source file
    :param weight: precalculated interpolation weight
    :param target: scrip grid describing target grid
    :param variables: list of variables to extract from the file, paired if the name in the file is different
    :param all_files: list of all files used as a source for these variables
    :param outfile: output filename, if none is generated as tempfile
    :param options: cdo options
    :param verbose: if events are printed
    :return: output filename
    """

    if len(lsm_file) > 0:
        # Read land-sea mask:
        nc = netCDF4.Dataset(lsm_file,'r')
        lsm = nc.variables[lsm_var][0,:]
        nc.close()
        # Define land-sea mask (True on the land):
        lsm_mask = (lsm>=lsm_limit)
    # First we deal with the variables which are interpolated using one
    # single input variable:
    varlist1 = []
    for x in variables:
        if isinstance(x['in'], list):
            varlist1 += x['in']
        else:
            varlist1.append(x['in'])
    if len(lsm_file) > 0:
        # Mask ERA variables: first we need to make sure that we are working
        # on a local file, not on one in /net/sea/work/datasets:
        file_out = f"{outdir}/{os.path.basename(file)}"
        if not os.path.exists(file_out):
            print(f"   copy {file} --> {file_out}")
            sys.stdout.flush()
            var_list = ",".join(varlist1)
            #cmd = f"{Global_settings.ncks} -v {var_list} {file} {file_out}"
            # Copy and convert ERA file to single precision:
            subprocess.run([Global_settings.cdo, "-b", "F32", f"select,name={var_list}", file, f"{outdir}/era_tmp.nc"], check=True)
            # Check if the ERA file contains scale_factor or add_offset attributes in the
            # variables we are interested in:
            nc = netCDF4.Dataset(f"{outdir}/era_tmp.nc",'r')
            scale_att_found = False
            for var in varlist1:
                vobj = nc.variables[var]
                if hasattr(vobj,"scale_factor") or hasattr(vobj,"add_offset"):
                    scale_att_found = True
            nc.close()
            if scale_att_found:
                # Unpack data to get rid of scale_factor and add_offset attributes:
                print('   unpack data using ncpdq')
                sys.stdout.flush()
                #subprocess.run([Global_settings.nccopy, "-k", "4", "-v", var_list, f"{outdir}/era_tmp1.nc", f"{outdir}/era_tmp2.nc"], check=True)
                subprocess.run([Global_settings.ncpdq, "--unpack", f"{outdir}/era_tmp.nc", file_out], check=True)
                subprocess.run(["/usr/bin/rm", "-f", f"{outdir}/era_tmp.nc"], check=True)
                # Make sure missing_value attributes are correct: ncpdq has adapted the _FillValue
                # attributes but not the missing_value
                nc = netCDF4.Dataset(file_out,'r+')
                for var in nc.variables:
                    vobj = nc.variables[var]
                    if hasattr(vobj,"missing_value"):
                        vobj.missing_value = getattr(vobj, "_FillValue")
                nc.close()
            else:
                subprocess.run(["/usr/bin/mv",f"{outdir}/era_tmp.nc", file_out], check=True)
        nc = netCDF4.Dataset(file_out,'a')
        #nt = len(nc.dimensions["time"])
        # Check for _FillValue and missing_value attributes which are not NaN:
        vlist_fillval = []
        for v in varlist1:
            vobj = nc.variables[v]
            if hasattr(vobj,'_FillValue'):
                fillval = getattr(vobj,"_FillValue")
                if not np.isnan(fillval):
                    vlist_fillval.append(v)
        if len(vlist_fillval) > 0:
            print("   set missing values to NaN")
            sys.stdout.flush()
            nc.close()
            # Replace missing values by NaN and adapt _FillValue atrributes:
            cmd_ncap2 = ""
            for v in vlist_fillval:
                cmd_ncap2 += f"{v}={v}; {v}.change_miss(nanf);"
            tmp_file = f"{outdir}/{str(uuid.uuid4())}"
            subprocess.run([Global_settings.ncap2, "-s", cmd_ncap2,file_out,tmp_file], check=True)
            subprocess.run(["/usr/bin/mv",tmp_file,file_out], check=True)
            # Reopen Netcdf file:
            nc = netCDF4.Dataset(file_out,'a')
            # Adapt missing_value attributes:
            for v in varlist1:
                vobj = nc.variables[v]
                if hasattr(vobj,'missing_value'):
                    vobj.missing_value = np.float32(np.nan)
        print(f"   apply land-sea mask from {lsm_file}")
        for v in varlist1:
            print(f"   apply ERA lsm to {v}")
            sys.stdout.flush()
            vobj = nc.variables[v]
            data = vobj[:]
            data[:, lsm_mask] = np.float32(np.nan)
            vobj[:] = data
            vobj.missing_value = np.float32(np.nan)
        nc.close()
        # Fill in missing values:
        print("   extrapolate to land")
        sys.stdout.flush()
        os.replace(file_out, f"{outdir}/tmp.nc")
        #cmd = f"{Global_settings.cdo} -s -P {Global_settings.cdo_omp_nthreads} -fillmiss2  {outdir}/tmp.nc {file_out}"
        #cmd = f"{Global_settings.cdo} -s -P 16 -fillmiss2  {outdir}/tmp.nc {file_out}"
        start_time = time.time()
        # sys.path.append(os.path.dirname(__file__)+"/../../../c++/build")
        # import roms_forcing_utils
        # roms_frc = roms_forcing_utils.ROMS_frc("", target, outdir, True)
        # roms_frc.era_extrapolate_land(f"{outdir}/tmp.nc", varlist1, Global_settings.cdo_omp_nthreads)
        subprocess.run([Global_settings.cdo, "-s", "fillmiss2", f"{outdir}/tmp.nc", file_out], check=True)
        print(f"      time for extrapolation [s]: {time.time()-start_time:.1f}")
        try:
            os.remove(f"{outdir}/tmp.nc")
        except OSError:
            pass
        print("   extrapolation done")
    else:
        # No masking of ERA data to be done: use whatever "file" points to:
        file_out = file
    #varlist1 = [x['in'] for x in variables if 'in' in x]
    varlist2 = [x['out'] for x in variables if 'expr' in x]
    renames = [x for x in variables if 'in' in x and x['in'] != x['out']]
    n_threads = Global_settings.cdo_omp_nthreads
    outfiles = []
    if outfile is None:
        #outfile = cdo.remap(target + ',' + weight, input=(' -selname,' + ','.join(varlist1) + ' ' + file_out),
        #                    options=options)
        outfile = f"{outdir}/{str(uuid.uuid4())}"
    if not os.path.exists(outfile):
        print("   remap forcing data to ROMS grid")
        start_time = time.time()
        cmd = f"{Global_settings.cdo} -s -P {n_threads} remap,{target},{weight} -selname,{','.join(varlist1)} {file_out} {outfile}"
        subprocess.run(cmd, shell=True, check=True)
        print(f"      time for horizontal regridding [s]: {time.time()-start_time:.1f}")
        outfiles.append(outfile)

    # Now deal with those variables which require a cdo expression rather than just one
    # single input variable:
    # the cdo expression is assumed to be x['expr'] if x in variables:
    # NOT necessary: is done by adjustments
    # for x in variables:
    #     if not 'expr' in x:
    #         continue
    #     outfile = f"{outdir}/{str(uuid.uuid4())}"
    #     cmd = f"{Global_settings.cdo} -s -P {n_threads} remap,{target},{weight} -expr,'{x['expr']}' {file_out} {outfile}"
    #     print(cmd)
    #     subprocess.run(cmd, shell=True, check=True)
    #     outfiles.append(outfile)

    # Merge datasets into one file if necessary:
    if len(outfiles) > 1:
        all_equal = True
        for i in range(1,len(outfiles)):
            if outfiles[0] != outfiles[i]:
                all_equal = False
        if all_equal:
            outfile = outfiles[0]
        else:
            #outfile = cdo.merge(input=(' '.join(outfiles)), options=options)
            outfile = f"{outdir}/{str(uuid.uuid4())}"
            subprocess.run([Global_settings.cdo,"-s","merge",' '.join(outfiles),outfile], check=True)
    elif len(outfiles) == 1:
        outfile = outfiles[0]

    # Rename output variable if necessary, and set some varible attributes:
    method = os.path.split(weight)[1].split("_")[0]
    with netCDF4.Dataset(outfile, mode='r+') as nc_file:
        # dims = nc_file.variables[varlist[0]].dimensions
        # if dims[-1] != "xi_rho":
        #     nc_file.renameDimension(dims[-1], "xi_rho")
        # if dims[-2] != "eta_rho":
        #     nc_file.renameDimension(dims[-2], "eta_rho")
        # if "time" not in dims[-3] and dims[-3] != "depth":
        #     nc_file.renameDimension(dims[-3], "depth")
        for name in varlist1+varlist2:
            if 'name' in nc_file.variables:
                v: netCDF4.Variable = nc_file.variables[name]
                v.setncattr('files', all_files)
                v.setncattr('h_interp_mtd', method)
        # for var in renames:
        #     if var['in'] in nc_file.variables:
        #         #if verbose:
        #         #    print('File '+outfile+': renaming ' + var['in'] + ' to ' + var['out'])
        #         nc_file.renameVariable(oldname=var['in'], newname=var['out'])
    #if verbose > 0:
    #    print('Sources and renaming complete!')
    return outfile


def calculate_weights(cdo, target_dir: str, group: dict, group_index: int, roms_grid: str, options: str,
                      verbose: int, in_file: str = ''):
    """
    calculates weights to use to interpolate
    :param cdo: cdo object
    :param target_dir: location to store weights under
    :param sources: list of dictionaries describing sources
    :param roms_grid: grid file of target grid
    :param options: cdo options
    :param verbose: whether a string is printed when this function is called
    :return: None
    """
    # Unfortunately, calling 'cdo genbil ...' twice with the same cdo object gives an error.
    # So generate a new cdo object:
    del cdo
    import cdo as mcdo
    cdo = mcdo.Cdo(debug=(verbose>2))
    cdo_mk_weight_mtd = {
        'bil': cdo.genbil,
        'bic': cdo.genbic,
        'nn': cdo.gennn,
        'dis': cdo.gendis,
        'con': cdo.gencon,
        'laf': cdo.genlaf
    }
    mtd_name = group['interpolation_method']
    # Path of file containing the weights:
    weight_name = os.path.join(target_dir, mtd_name + '_weight_g' + str(group_index) + ".nc")
    group['weight'] = weight_name
    if os.path.exists(weight_name):
        return
    if verbose > 0:
        print("Calculating horiz. interpolation weights and storing them in: " + weight_name)
    mtd = cdo_mk_weight_mtd[mtd_name]
    if 'files' in group:
        files = group['files']
        a_file = files if isinstance(files, str) else files[0]
    else:
        a_file = in_file
    mtd(roms_grid, input=a_file, output=weight_name, options=options)


"""
This function was originally intended for interpolating grb files. 
Due to time constraints, this functionality was not finished.
Please fix this function to have the same functionaity as above.
Please view romspy/data_adjustment/grb_tools.py 
"""
#
# def __interpolate_grb(file: str, variables: list, weight: str, all_files: str, outfile: str):
#     """
#     Convert a .grd file to a .nc file using cdo while extracting the variable and interpolating
#     :param file: file to convert
#     :param variables: list of variables to extract from the file, paired if the name in the file is different
#     :param weight: precalculated interpolation weight
#     :param all_files: list of all files used as a source for these variables
#     :param outfile: output filename
#     """
#
#     # If era numbering changes change the parameter table
#     my_table: str = grb_tools.get_para_table_file_vars(grb_tools.para_table_era_interim, variables, file, all_files)
#     varlist = [a['name'] for a in grb_tools.para_table_era_interim if a['out_name'] in variables]
#     # -selname select the variable by name from file first
#     # -setpartabn rename the variable and set attributes according to table second
#     # remap weight is the argument, no need to specify method
#     # -f file type
#     # -P use 8 cores
#     ROMS.my_cdo.remap(self.scrip_grid + ',' + weight,
#                       input='-setpartabn,' + my_table + ' -selname,' + ','.join(varlist) + ' ' + file,
#                       options=self.options, output=outfile)
#     os.remove(my_table)
