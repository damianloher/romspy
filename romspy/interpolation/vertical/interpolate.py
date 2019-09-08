import netCDF4
import numpy as np
import os

"""
Author: Nicolas Munnich
License: GNU GPL2+
"""


def vert_interpolate(cdo, gen_fun, apply_fun, weight_extra_len: int, file: str, outfile: str, weight_dir: str,
                     variables: list, z_levels: np.ndarray, source: dict, options: str,
                     verbose: bool = False) -> str:
    """
    Interpolate a number of 4D variables in a file on the third to last axis.
    For example, a variable with the dimensions (time, depth, eta, xi) will be interpolated on the depth axis
    for each point in time.
    :param cdo: cdo object
    :param gen_fun: C external function which generates weights
    :param apply_fun: C external function which applies weights
    :param weight_extra_len: length of the final dim in weights
    :param file: File which contains the variables
    :param outfile: Name of output file
    :param weight_dir: directory to store weights in
    :param variables: Variables to be vertically interpolated. Must all be on the same grid.
    :param z_levels: 3D array of depths at each point
    :param vertical_weights: dictionary of vertical weights and their creation circumstances
    :param options: cdo options
    :param verbose: if runtime info should be printed to consoles
    :return: File with variables vertically interpolated
    """

    # weight_file = source.get('vertical_weight', None)

    # Produce weight if it does not already exist or replace weight if ssh is found
    # if weight_file is None:
    weight_file = gen_vert_weight(gen_fun, weight_extra_len, weight_dir, file, variables, z_levels, verbose)
    # source['vertical_weight'] = weight_file

    # Interpolate with the corresponding weight
    return apply_vert_weights(cdo, apply_fun, weight_file, file, outfile, variables, options, verbose)


def gen_vert_weight(gen_fun, weight_extra_len: int, weight_dir: str, file: str, variables: list, z_levels: np.ndarray,
                    verbose: bool = False) -> str:
    """
   Generate weights for vertical interpolation routines.
    :param gen_fun: Vertical interpolation routine written in C. Must have same arguments as create_weights
                    found in romspy/interpolation/vertical/linear
    :param weight_extra_len: length of the extra dimension added to weights
    :param weight_dir: Directory to store weights in
    :param file: input file to create weights for
    :param variables: name of the variables used by this weight, used to create weight name
    :param z_levels: 3D array of depth points
    :param verbose: whether to print whenever this function is called
    :return: name of weight file
    """
    if verbose:
        print("Making weight")
    with netCDF4.Dataset(file, mode='r') as my_file:
        # Get Z to interpolate to and its length
        z_levels = np.require(z_levels, requirements=["C", "O", "A"], dtype=np.float32)
        z_slice = np.prod(z_levels.shape[1:], dtype=np.ulonglong)
        z_len = np.ulonglong(z_levels.shape[0])
        # Get weight array which will store results
        weight_arr = np.require(np.empty(tuple(z_levels.shape + (weight_extra_len,)), dtype=np.float32),
                                requirements=['C', "W", "O", "A"])  # 2 here is [index, index percentage] see c code

        # Get depth and depth length
        depth = np.require(my_file.variables["depth"][:], dtype=np.float32, requirements=["C", "O", "A"])
        depth_len = np.ulonglong(depth.size)
        if verbose:
            print("Calling external C routine for generating weights.")
        # GET WEIGHTS
        gen_fun(depth, depth_len, z_levels, z_slice, z_len, weight_arr)
        if verbose:
            print("C routine has been completed.")
        # Get name to save weights under
        prefix = "vertical_" + "_".join(variables) + "_"
        path = os.path.split(file)[1]
        weight_file_name = os.path.join(weight_dir, prefix + path + ".npy")
        # print(weight_arr)
        # Save weights
        np.save(weight_file_name, weight_arr)
    return weight_file_name


def apply_vert_weights(cdo, apply_fun, weight_file: str, file: str, outfile: str, variables: list, options: str,
                       verbose: bool):
    if verbose:
        print("Interpolating vertically with weights: " + " ".join(variables))
    # Get weight
    weight_data = np.load(weight_file)
    # Get a temporary filepath
    split = os.path.split(file)
    temp_out_path = os.path.join(split[0], "temp_vert_" + split[1])
    # Open input file
    with netCDF4.Dataset(file, mode="r+") as in_file:
        # Get the dimension names used in the variables
        dims = in_file.variables[variables[0]].dimensions
        var_dims: list = [x if x != "depth" else "s_rho" for x in dims]
        # Get dimension lengths
        var_dim_lens = weight_data.shape[:-1]
        time_len = len(in_file.dimensions["time"])
        # Create destination array
        var_out_arr = np.require(np.zeros(tuple(var_dim_lens), dtype=np.float32), requirements=["C", "W", "O", "A"])
        arr_len = np.ulonglong(var_out_arr.size)

        with netCDF4.Dataset(temp_out_path, mode="w") as dest_file:
            # Create dimensions
            dest_file.createDimension("time", time_len)
            dest_file.createDimension("s_rho", var_dim_lens[0])
            dest_file.createDimension(var_dims[-2], var_dim_lens[1])
            dest_file.createDimension(var_dims[-1], var_dim_lens[2])
            for var in variables:
                # create variable and set attributes
                new_var: netCDF4.Variable = dest_file.createVariable(var, 'f', tuple(var_dims))
                var_obj = in_file.variables[var]
                var_attrs = {x: str(var_obj.getncattr(x)) for x in var_obj.ncattrs() if x != "_FillValue"}
                new_var.setncatts(var_attrs)
                # Time step wise interpolation
                if verbose:
                    print("Interpolating " + var + " vertically with external C routine timestep wise.")
                for i in range(time_len):
                    var_data = np.ascontiguousarray(var_obj[i], dtype=np.float32)
                    apply_fun(weight_data, var_data, var_out_arr, arr_len)
                    new_var[i] = np.require(var_out_arr, dtype=np.float32)
                if verbose:
                    print("Finished interpolating " + var + " vertically with external C routine.")
                in_file.renameVariable(var, "tmp_" + var)
    # Merge into dataset
    # doesn not work CDO limitations! cdo.replace(input=file + " " + temp_out_path, output=outfile, options=options)
    temp_merge = cdo.merge(input=file + " " + temp_out_path, options=options)
    if outfile is not None:
        cdo.delname(",".join(["tmp_" + var for var in variables]), input=temp_merge, output=outfile, options=options)
    else:
        outfile = cdo.delname(",".join(["tmp_" + var for var in variables]), input=temp_merge, options=options)
    os.remove(temp_out_path)
    # dest_dir = os.path.split(outfile)
    # temp_file = os.path.join(dest_dir[0], "temp_file.nc")
    # cdo.copy(input=temp_out_path, output = temp_file, options = options)
    return outfile
