from itertools import count
import netCDF4
import os

"""
Author: Nicolas Munnich
License: GNU GPL2+
"""


def cdo_interpolate(cdo, file: str, weight: str, target_grid: str, variables: list, all_files: str, options: str,
                    outfile: str = None, verbose=False) -> str:
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
    if verbose:
        print('Horizontally interpolating ' + " ".join([x['out'] for x in variables]) + " from " + file)

    import netCDF4
    try:
        nc = netCDF4.Dataset(file,'r')
        nc.close()
        return __interpolate_nc(cdo, file, weight, target_grid, variables, all_files, outfile, options, verbose)
    except OSError:
        #extension = "." + file[::-1].split(".")[0][::-1]  # get file extension
        #raise ValueError(extension + " source files are currently not supported. "
        raise ValueError("not a Netcdf file: " + file + "\n" +
                         "If you wish to add this functionality, " +
                         "please look under romspy/interpolation/horizontal/interpolate.py")
    # Room for additional/alternative conversion functions if they appear necessary


def __interpolate_nc(cdo, file: str, weight: str, target: str, variables: list, all_files: str, outfile: str,
                     options: str, verbose: bool) -> str:
    """

    interpolated and extract variable
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

    # First we deal with the variables which are interpolated using one
    # single input variable:
    varlist1 = [x['in'] for x in variables if 'in' in x]
    varlist2 = [x['out'] for x in variables if 'expr' in x]
    renames = [x for x in variables if 'in' in x and x['in'] != x['out']]
    # -selname select the variable by name from file first
    # remap weight is the argument, no need to specify method
    # -f file type
    # -P use 8 cores
    outfiles = []
    if outfile is None:
        outfile = cdo.remap(target + ',' + weight, input=(' -selname,' + ','.join(varlist1) + ' ' + file),
                            options=options)
    else:
        cdo.remap(target + ',' + weight, input=(' -selname,' + ','.join(varlist1) + ' ' + file),
                  options=options, output=outfile)
    outfiles.append(outfile)

    # Now deal with those variables which require a cdo expression rather than just one
    # single input variable:
    # the cdo expression is assumed to be x['expr'] if x in variables:
    for x in variables:
        if not 'expr' in x:
            continue
        outfile = cdo.remap(target + ',' + weight, input=(" -expr,'{}'".format(x['expr']) + ' ' + file),
                            options=options)
        outfiles.append(outfile)

    # Merge datasets into one file if necessary:
    if len(outfiles) > 1:
        outfile = cdo.merge(input=(' '.join(outfiles)), options=options)

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
        for var in renames:
            if var['in'] in nc_file.variables:
                #if verbose:
                #    print('File '+outfile+': renaming ' + var['in'] + ' to ' + var['out'])
                nc_file.renameVariable(oldname=var['in'], newname=var['out'])
    if verbose:
        print('Sources and renaming complete!')
    return outfile


def calculate_weights(cdo, target_dir: str, sources: list, scrip_grid: str, options: str,
                      verbose: bool, in_file: str = ''):
    """
    calculates weights to use to interpolate
    :param cdo: cdo object
    :param target_dir: location to store weights under
    :param sources: list of dictionaries describing sources
    :param scrip_grid: grid to interpolate onto
    :param options: cdo options
    :param verbose: whether a string is printed when this function is called
    :return: None
    """
    # Unfortunately, calling 'cdo genbil ...' twice with the same cdo object gives an error.
    # So generate a new cdo object:
    del cdo
    import cdo as mcdo
    # cdo_mk_weight_mtd = {
    #     'bil': cdo.genbil,
    #     'bic': cdo.genbic,
    #     'nn': cdo.gennn,
    #     'dis': cdo.gendis,
    #     'con': cdo.gencon,
    #     'con2': cdo.gencon2,
    #     'laf': cdo.genlaf
    # }
    if verbose:
        print("Calculating weights. Target directory: " + target_dir)
    for group, index in zip(sources, count()):
        cdo = mcdo.Cdo(debug=verbose)
        #print('new cdo object generated!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!')
        cdo_mk_weight_mtd = {
            'bil': cdo.genbil,
            'bic': cdo.genbic,
            'nn': cdo.gennn,
            'dis': cdo.gendis,
            'con': cdo.gencon,
            'con2': cdo.gencon2,
            'laf': cdo.genlaf
        }
        mtd_name = group['interpolation_method']
        weight_name = os.path.join(target_dir, mtd_name + '_weight_' + str(index) + ".nc")
        mtd = cdo_mk_weight_mtd[mtd_name]
        if 'files' in group:
            files = group['files']
            a_file = files if isinstance(files, str) else files[0]
        else:
            a_file = in_file
        mtd(scrip_grid, input=a_file, output=weight_name, options=options)
        group['weight'] = weight_name


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
