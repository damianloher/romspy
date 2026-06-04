import os
import netCDF4
import numpy as np

"""
Author: Nicolas Munnich
License: GNU GPL2+
"""


def adjust_vectors(cdo, in_file, target_grid, variables, options, verbose=True, out_file=None) -> str:
    """
    Shifts and rotates variable pairs
    :param cdo
    :param in_file:
    :param target_grid:
    :param variables:
    :param options:
    :param verbose:
    :param out_file:
    :return:
    """

    split = os.path.split(in_file)
    temp_out_path = os.path.join(split[0], "temp_" + split[1])
    with netCDF4.Dataset(target_grid, mode='r') as target:
        angle = target.variables["angle"]
        if angle.units == "degrees":
            angle = np.deg2rad(angle[:])
        else:
            angle = angle[:]
    with netCDF4.Dataset(in_file, mode='r+') as _in:
        with netCDF4.Dataset(temp_out_path, mode="w") as _out:
            if not "xi_u" in _out.dimensions:
                if "xi_rho" in _in.dimensions:
                    _out.createDimension("xi_u", len(_in.dimensions["xi_rho"]) - 1)
                elif "xi_u" in _in.dimensions:
                    _out.createDimension("xi_u", len(_in.dimensions["xi_u"]))
                else:
                    _out.createDimension("xi_u", len(_in.dimensions["xi_v"] - 1))
            if not "eta_u" in _out.dimensions:
                if "eta_rho" in _in.dimensions:
                    _out.createDimension("eta_u", len(_in.dimensions["eta_rho"]))
                elif "eta_u" in _in.dimensions:
                    _out.createDimension("eta_u", len(_in.dimensions["eta_u"]))
                else:
                    _out.createDimension("eta_u", len(_in.dimensions["eta_v"]) + 1)
            if not "xi_v" in _out.dimensions:
                if "xi_rho" in _in.dimensions:
                    _out.createDimension("xi_v", len(_in.dimensions["xi_rho"]))
                elif "xi_u" in _in.dimensions:
                    _out.createDimension("xi_v", len(_in.dimensions["xi_u"]) + 1)
                else:
                    _out.createDimension("xi_v", len(_in.dimensions["xi_v"]))
            if not "eta_v" in _out.dimensions:
                if "eta_rho" in _in.dimensions:
                    _out.createDimension("eta_v", len(_in.dimensions["eta_rho"]) - 1)
                elif "eta_u" in _in.dimensions:
                    _out.createDimension("eta_v", len(_in.dimensions["eta_u"]) - 1)
                else:
                    _out.createDimension("eta_v", len(_in.dimensions["eta_v"]))
            for dim in _in.variables[variables[0][0]].dimensions:
                dim_len = len(_in.dimensions[dim])
                if not dim in _out.dimensions:
                    _out.createDimension(dim, dim_len)
                if dim == "depth":
                    d_obj = _out.createVariable("depth", 'd', ("depth",))
                    d_obj[:] = _in.variables["depth"][:]
                    d_obj.setncattr("units", "meters")
                    d_obj.setncattr("positive", "down")
            for u, v in variables:
                if verbose:
                    print("Making vectors: (" + u + "," + v + ")")
                u_obj, v_obj = _in.variables[u], _in.variables[v]
                dims = list(u_obj.dimensions)
                u_dims, v_dims = dims.copy(), dims
                u_dims[-1] = "xi_u"
                u_dims[-2] = "eta_u"
                v_dims[-1] = "xi_v"
                v_dims[-2] = "eta_v"
                is_3d = len(u_dims) > 3

                time_length = len(_in.dimensions[dims[0]])
                new_u: netCDF4.Variable = _out.createVariable(u, 'f', tuple(u_dims))
                new_v: netCDF4.Variable = _out.createVariable(v, 'f', tuple(v_dims))
                new_u.setncatts({x: u_obj.getncattr(x) for x in u_obj.ncattrs()})
                new_v.setncatts({x: v_obj.getncattr(x) for x in v_obj.ncattrs()})
                for t in range(time_length):
                    u_contents, v_contents = u_obj[t], v_obj[t]
                    # Rotate u and v components on rho grid:
                    cosa = np.cos(angle)
                    sina = np.sin(angle)
                    u_turned = u_contents * cosa + v_contents * sina
                    v_turned = v_contents * cosa - u_contents * sina
                    # Interpolate to u and v grid, respectively:
                    if is_3d:
                        u_contents = 0.5 * (u_turned[:, :, 1:] + u_turned[:, :, :-1])
                        v_contents = 0.5 * (v_turned[:, 1:, :] + v_turned[:, :-1, :])
                    else:
                        u_contents = 0.5 * (u_turned[:, 1:] + u_turned[:, :-1])
                        v_contents = 0.5 * (v_turned[1:, :] + v_turned[:-1, :])

                    new_u[t] = u_contents
                    new_v[t] = v_contents
                _in.renameVariable(u, "tmp_" + u)
                _in.renameVariable(v, "tmp_" + v)

    #t_name = cdo.merge(input=in_file + " " + temp_out_path, options=options)
    #if not os.path.exists(t_name):
    if in_file != temp_out_path:
        # Execute cdo command directly:
        t_name = temp_out_path + "_2"
        cmd = f"/usr/local/bin/cdo {options} -merge {in_file} {temp_out_path} {t_name}"
        if verbose:
            print(cmd)
        os.system(cmd)
    else:
        t_name = in_file
    if out_file is not None:
        cdo.delname(",".join(["tmp_" + u + ",tmp_" + v for u, v in variables]), input=t_name, output=out_file,
                    options=options)
    else:
        out_file = cdo.delname(",".join(["tmp_" + u + ",tmp_" + v for u, v in variables]), input=t_name,
                               options=options)
    os.remove(temp_out_path)
    return out_file


def shift(h: np.ndarray, grid_type: int):
    """
    shift a 2d grid
    :param h:
    :param grid_type: 0 if rho-rho, 1 if rho-u, 2 if v-rho
    :return:
    """
    if grid_type == 0:
        return h
    elif grid_type == 1:
        return 0.5 * (h[:, 1:] + h[:, :-1])
    elif grid_type == 2:
        return 0.5 * (h[1:, :] + h[:-1, :])
