from romspy.interpolation.vertical.levels import z_levels, sigma_stretch_cs, sigma_stretch_sc
from romspy.interpolation.shift_grid import shift
import netCDF4
import numpy as np

"""
Author: Nicolas Munnich
License: GNU GPL2+
"""


def uv_bar_adjustment(file: str, group_files: str, target_grid, verbose, h, layers, theta_s, theta_b, sigma_type, hc,
                      zeta, obc,
                      **kwargs):
    if verbose:
        print("Getting ubar and vbar")
    h = h
    sc = sigma_stretch_sc(layers, False)
    cs = sigma_stretch_cs(theta_s, theta_b, sc, sigma_type)
    hc = hc
    zeta = zeta
    delta_z_levs_u = z_levels(shift(h, 1), sc, cs, hc, shift(zeta, 1), sigma_type=3, verbose=verbose)
    delta_z_levs_u = delta_z_levs_u[1:] - delta_z_levs_u[:-1]
    delta_z_levs_u_inv = 1 / sum(delta_z_levs_u)
    delta_z_levs_v = z_levels(shift(h, 2), sc, cs, hc, shift(zeta, 2), sigma_type=3, verbose=verbose)
    delta_z_levs_v = delta_z_levs_v[1:] - delta_z_levs_v[:-1]
    delta_z_levs_v_inv = 1 / sum(delta_z_levs_v)
    my_file: netCDF4.Dataset
    with netCDF4.Dataset(file, mode='r+') as my_file:
        time_len = len(my_file.dimensions["time"])
        u: netCDF4.Variable = my_file.variables["u"]
        v: netCDF4.Variable = my_file.variables["v"]
        ubar_dims = (u.dimensions[0], u.dimensions[2], u.dimensions[3])
        vbar_dims = (v.dimensions[0], v.dimensions[2], v.dimensions[3])
        ubar: netCDF4.Variable = my_file.createVariable('ubar', 'f', ubar_dims)
        vbar: netCDF4.Variable = my_file.createVariable('vbar', 'f', vbar_dims)
        u.setncattr('long_name', "u-velocity component")
        v.setncattr('long_name', "v-velocity component")
        ubar_attrs = {x: u.getncattr(x) for x in u.ncattrs()}
        ubar_attrs['long_name'] = "vertically integrated u-velocity component"
        vbar_attrs = {x: v.getncattr(x) for x in v.ncattrs()}
        vbar_attrs['long_name'] = "vertically integrated v-velocity component"
        ubar.setncatts(ubar_attrs)
        vbar.setncatts(vbar_attrs)

        with netCDF4.Dataset(target_grid) as grd:
            rmask = grd.variables['mask_rho'][:]
            pn = grd.variables['pn'][:]
            pm = grd.variables['pm'][:]
        for t in range(time_len):
            ubar_vals = np.sum(np.multiply(u[t], delta_z_levs_u), 0)
            vbar_vals = np.sum(np.multiply(v[t], delta_z_levs_v), 0)
            ubar_vals, vbar_vals = get_obcvolcons(ubar_vals, vbar_vals, pm, pn, rmask, obc, verbose)
            ubar_vals *= delta_z_levs_u_inv
            vbar_vals *= delta_z_levs_v_inv
            ubar[t] = np.where(u[t] > -9.9e32, ubar_vals, u[t])
            vbar[t] = np.where(v[t] > -9.9e32, vbar_vals, v[t])


def get_obcvolcons(ubar, vbar, pm, pn, rmask, obc, verbose):
    """
    Enforce integral flux conservation around the domain
    :param ubar: u averaged over depth
    :param vbar: v averaged over depth
    :param pm: curvilinear coordinate metric in xi
    :param pn: curvilinear coordinate metric in eta
    :param rmask: rho-grid mask
    :param obc: open boundary conditions, (1=open 0=closed, [S E N W])
    :param verbose: whether to print runtime information
    :return: 
    """
    umask = 2 * shift(rmask, 1)
    vmask = 2 * shift(rmask, 2)
    dy_u = 2 * umask / (2 * shift(pn, 1))
    dx_v = 2 * vmask / (2 * shift(pm, 2))
    udy = ubar * dy_u
    vdx = vbar * dx_v
    flux = obc[0] * np.nansum(vdx[1:, 0]) - obc[1] * np.nansum(udy[-1, 1:]) - obc[2] * np.nansum(vdx[1:, -1]) + \
           obc[3] * np.nansum(udy[0, 1:])
    cross = obc[0] * np.nansum(dx_v[1:, 0]) + obc[1] * np.nansum(dy_u[-1, 1:]) + obc[2] * np.nansum(dx_v[1:, -1]) + \
            obc[3] * np.nansum(dy_u[0, 1:])
    vcorr = flux / cross
    if verbose:
        print("Flux correction: " + str(vcorr))
    vbar[:, 0] = obc[0] * (vbar[:, 0] - vcorr)
    ubar[-1, :] = obc[1] * (ubar[-1, :] + vcorr)
    vbar[:, -1] = obc[2] * (vbar[:, -1] + vcorr)
    ubar[0, :] = obc[3] * (ubar[0, :] - vcorr)

    ubar = ubar * umask
    vbar = vbar * vmask
    return ubar, vbar


clim_adjustments = [
    {'out_var_names': {'ubar', 'vbar'}, 'in_var_names': {'u', 'v'},
     'func': uv_bar_adjustment}
]
