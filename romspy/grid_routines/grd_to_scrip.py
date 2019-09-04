import netCDF4
import numpy as np
import xarray as xr
import math
import os

"""
Author: Nicolas Munnich
License: GNU GPL2+
"""


# TODO: Investigate bugged con
def from_nc(filepath: str) -> str:
    """
    Make a SCRIP grid file from a grid file
    :param filepath: grid file
    :return: SCRIP grid file
    """
    given_grid = netCDF4.Dataset(filepath)
    lon = given_grid.variables['lon_rho'][:]
    lat = given_grid.variables['lat_rho'][:]
    mask = given_grid.variables['mask_rho'][:]
    path = os.path.split(filepath)
    if not os.path.isdir(os.path.join(path[0], "scrip_grid")):
        os.mkdir(os.path.join(path[0], "scrip_grid"))
    out_string = os.path.join(path[0], "scrip_grid", 'scrip_' + path[1])
    title = filepath
    given_grid.close()

    return scrip_grd_maker(lon, lat, mask, out_string, title)


def scrip_grd_maker(lon: np.ndarray, lat: np.ndarray, mask: np.ndarray, outfile: str, title: str) -> str:
    """
    Make a SCRIP grid
    :param lon: 2d np array of degree longitudes of the grid_routines centers
    :param lat: 2d np array of degree latitudes of the grid_routines centers
    :param mask: 2d np array of land mask (1 ocean point, 0 land point)
    :param outfile: output NETCDF filename
    :param title: title string to add
    :return: outfile for consistency
    """
    assert lon.shape == lat.shape == mask.shape
    # grid_size = lon.shape[0] * lon.shape[1]
    # grid_corners = 4
    # grid_rank = 2

    deg2rad = math.pi / 180
    lon = lon * deg2rad
    lat = lat * deg2rad

    # Compute corner position
    # lon: half-way, lat half-way weighted by cos(lat) to account for smaller area at higher lat

    lon_xy = 0.5 * (fbar_x(fbar_y(lon)) + fbar_y(fbar_x(lon)))
    latbar_y = fbar_y(lat, np.cos(lat))
    latbar_x = fbar_x(lat, np.cos(lat))
    lat_xy = 0.5 * (fbar_x(latbar_y, np.cos(latbar_y)) + fbar_y(latbar_x, np.cos(latbar_x)))

    corner_lon_data = np.stack(
        (np.ravel(lon_xy[:-1, :-1], order='F'), np.ravel(lon_xy[1:, :-1], order='F'),
         np.ravel(lon_xy[1:, 1:], order='F'), np.ravel(lon_xy[:-1, 1:], order='F')))
    corner_lat_data = np.stack(
        (np.ravel(lat_xy[:-1, :-1], order='F'), np.ravel(lat_xy[1:, :-1], order='F'),
         np.ravel(lat_xy[1:, 1:], order='F'), np.ravel(lat_xy[:-1, 1:], order='F')))
    pi2 = math.pi / 2
    corner_lat_data = np.where(corner_lat_data < pi2, corner_lat_data, pi2)
    corner_lat_data = np.where(corner_lat_data > -pi2, corner_lat_data, -pi2)

    """
    Required format:
    global attributes:
    title = SCRIP grid_routines input file for: <inputfile>

    dims:
    grid_size
    grid_corners
    grid_rank

    variables:
    grid_dims (grid_rank) = lon.shape
    grid_center_lon (grid_size) = rad_lon
        @units = radians
    grid_center_lat (grid_size) = rad_lat
        @units = radians
    grid_imask (grid_size) = mask
    grid_corner_lon (grid_corners, grid_size) = corner_lon_data
        @units = radians
    grid_corner_lat (grid_corners, grid_size) = corner_lat_data
        @units = radians
    """
    grid_dims = xr.DataArray(np.array([lon.shape[1], lon.shape[0]]), dims=['grid_rank'], name='grid_dims')
    grid_center_lon = xr.DataArray(lon.flatten(), dims=['grid_size'], attrs={'units': 'radians'},
                                   name='grid_center_lon')
    grid_center_lat = xr.DataArray(lat.flatten(), dims=['grid_size'], attrs={'units': 'radians'},
                                   name='grid_center_lat')
    grid_imask = xr.DataArray(mask.astype('i').flatten(), dims=['grid_size'], name='grid_imask')
    grid_corner_lon = xr.DataArray(corner_lon_data.transpose(), dims=['grid_size', 'grid_corners'],
                                   name='grid_corner_lon', attrs={'units': 'radians'})
    grid_corner_lat = xr.DataArray(corner_lat_data.transpose(), dims=['grid_size', 'grid_corners'],
                                   name='grid_corner_lat', attrs={'units': 'radians'})
    scrip_grid = xr.Dataset(
        {'grid_dims': grid_dims, 'grid_center_lon': grid_center_lon, 'grid_center_lat': grid_center_lat,
         'grid_imask': grid_imask, 'grid_corner_lon': grid_corner_lon, 'grid_corner_lat': grid_corner_lat},
        attrs={'title': 'SCRIP input file for ' + title})

    scrip_grid.to_netcdf(outfile)
    return outfile


def fbar_x(f: np.ndarray, w: np.ndarray = None):
    """
    computes (weighted) neighbor mean in y direction and extrapolates half-way to the sides
    :param f: 2D numpy array to compute means of
    :param w: 2D numpy array or None, if no weights are present
    :return: weighted neighbor mean, shape (L, M) -> (L, M+1)
    """
    # Shift for offset
    front = np.pad(f, ((1, 0), (0, 0)), 'constant')
    back = np.pad(f, ((0, 1), (0, 0)), 'constant')
    if w is None:
        # adjust edges for extrapolation
        # 0.5 z + 0.5 x = 1.5 x - 0.5 y
        # z = 2x - y
        front[0] = 2 * front[1] - front[2]
        back[-1] = 2 * back[-2] - back[-3]
        return (front + back) * 0.5

    wbi = np.pad(1 / (w[:-1] + w[1:]), ((1, 1), (0, 0)), 'edge')
    # adjust edges for extrapolation
    # 2 * x * w_x + x * w_y - y * w_y = z * w_z + x * w_x
    # z * w_z = x * w_x + (x - y) * w_y
    # then set w_z to 1
    w_front = np.pad(w, ((1, 0), (0, 0)), 'constant', constant_values=1)
    w_back = np.pad(w, ((0, 1), (0, 0)), 'constant', constant_values=1)
    back[-1] = (back[-2] * w_back[-2] + (back[-2] - back[-3]) * w_back[-3])
    front[0] = (front[1] * w_front[1] + (front[1] - front[2]) * w_front[2])

    return (front * w_front + back * w_back) * wbi


def fbar_y(f: np.ndarray, w: np.ndarray = None):
    """
    computes (weighted) neighbor mean in y direction and extrapolates half-way to the sides
    :param f: 2D numpy array to compute means of
    :param w: 2D numpy array or None, if no weights are present
    :return: weighted neighbor mean, shape (L, M) -> (L, M+1)
    """
    # Shift for offset
    front = np.pad(f, ((0, 0), (1, 0)), 'constant')
    back = np.pad(f, ((0, 0), (0, 1)), 'constant')
    if w is None:
        # adjust edges for extrapolation
        # 0.5 z + 0.5 x = 1.5 x - 0.5 y
        # z = 2x - y
        front[:, 0] = 2 * front[:, 1] - front[:, 2]
        back[:, -1] = 2 * back[:, -2] - back[:, -3]
        return (front + back) * 0.5

    wbi = np.pad(1 / (w[:, :-1] + w[:, 1:]), ((0, 0), (1, 1)), 'edge')
    # adjust edges for extrapolation
    # 2 * x * w_x + x * w_y - y * w_y = z * w_z + x * w_x
    # z * w_z = x * w_x + (x - y) * w_y
    # then set w_z to 1
    w_front = np.pad(w, ((0, 0), (1, 0)), 'constant', constant_values=1)
    w_back = np.pad(w, ((0, 0), (0, 1)), 'constant', constant_values=1)
    back[:, -1] = (back[:, -2] * w_back[:, -2] + (back[:, -2] - back[:, -3]) * w_back[:, -3])
    front[:, 0] = (front[:, 1] * w_front[:, 1] + (front[:, 1] - front[:, 2]) * w_front[:, 2])

    return (front * w_front + back * w_back) * wbi
