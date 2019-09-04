"""
Main ROMS_TOOLS class.
@Author: Nicolas Munnich

Some choices were made which slightly reduce efficiency of code in name of readability.
These choices are mostly present when combining cdo and netCDF4 libraries.
A major limitation was produced by the file format of netCDF4 making deleting variables or adjusting their sizes
not possible to do so directly, but requires a roundabout approach to be used.
"""

import cdo
import netCDF4
import os
import numpy as np
from itertools import count

from romspy.verification import test_cdo, has_vertical, verify_sources
from .grid_routines import scrip_grid_from_nc
from .interpolation.vertical import sigma_stretch_sc, sigma_stretch_cs, get_z_levels
from .interpolation import Interpolator
from .settings import Settings

"""
Author: Nicolas Munnich
License: GNU GPL2+
"""


class PreProcessor:
    def __init__(self, target_grid: str, outfile: str, sources: list, setting: Settings,
                 vertical_options: dict = None, scrip_grid: str = None, vertical_weights: dict = None,
                 time_underscored: bool = False, processes: int = 8, file_type: str = 'nc4c',
                 keep_weights: bool = False, keep_z_clim: bool = False, verbose: bool = False):
        """
        Contains universal methods

        :param target_grid: grid_routines to be interpolated onto
        :param sources: Information about which variables should be taken from which files using
                which interpolation method onto which type of grid_routines and whether to interpolate them vertically
        :param setting: Setting object containing functions which should be
                        run on some variables to produce other variables
        :param scrip_grid: Scrip file on rho grid_routines to interpolate with. Will be created if not provided
        :param time_underscored: Whether the variable 'time' should be renamed to '#_time' for backwards compatibility
        :param processes: How many processes can be run in parallel
        :param file_type: output file type
        :type file_type:   'nc1' | 'nc2' | 'nc4' | 'nc4c'
                           recommended and default: 'nc4c'
        :param keep_weights: whether calculated weights should be kept
        :param verbose: whether text should be printed as the program is running
        """
        # cdo
        self.cdo = cdo.Cdo()
        test_cdo(self.cdo)

        # Sources
        verify_sources(sources, verbose)
        self.sources = sources
        self.target_grid = target_grid

        # Vertical interpolation information
        # Replace anything not passed in with default values
        vertical_options = {} if vertical_options is None else vertical_options
        vertical_options = {**{"theta_s": 7.0, "theta_b": 0.0, "layers": 32, "hc": 150.0, "tcline": 150.0,
                               "sigma_type": 3}, **vertical_options}
        self.vertical_options = vertical_options
        self.sc = sigma_stretch_sc(vertical_options["layers"], True)
        self.cs = sigma_stretch_cs(vertical_options["theta_s"],
                                   vertical_options["theta_b"], self.sc, vertical_options["sigma_type"])
        self.hc = vertical_options["hc"]
        self.vertical_information = self.get_vertical_information(vertical_options)

        # CDO options
        self.options = " -b F32 -f " + file_type + " -P " + str(processes)
        # Other Options
        self.verbose = verbose
        self.time_underscored = time_underscored
        self.outfile = outfile

        # adjustments which need to be made
        self.var_settings = setting
        setting.preprocessor = self

        # Get z_levels
        with netCDF4.Dataset(target_grid, mode='r') as my_grid:
            self.h = my_grid.variables['h'][:]
        self.zeta = np.zeros_like(self.h)
        self.all_z_levels = get_z_levels(self.h, self.sc, self.cs, self.hc, self.zeta)
        # To make zeta non-flat, call cdo_interpolate on the file holding ssh and then extract ssh like self.h

        # Interpolator
        self.scrip_grid = scrip_grid_from_nc(target_grid) if scrip_grid is None else scrip_grid
        self.interpolator = Interpolator(self.cdo, os.path.split(outfile)[0], sources, target_grid, self.scrip_grid,
                                         self.all_z_levels, self.options, vertical_weights, keep_weights, keep_z_clim,
                                         verbose)

        if verbose:
            print("Finished setup")

    def get_vertical_information(self, vertical_options: dict):
        return [
            {'name': 'theta_s', 'long_name': 'S-coordinate surface control parameter', 'datatype': 'f',
             'dimensions': 'one', 'units': '-', 'data': vertical_options["theta_s"]},
            {'name': 'theta_b', 'long_name': 'S-coordinate bottom control parameter', 'datatype': 'f',
             'dimensions': 'one', 'units': '-', 'data': vertical_options["theta_b"]},
            {'name': 'Tcline', 'long_name': 'S-coordinate surface/bottom layer width', 'datatype': 'f',
             'dimensions': 'one', 'units': 'meter', 'data': vertical_options["tcline"]},
            {'name': 'hc', 'long_name': 'S-coordinate critical depth', 'datatype': 'f',
             'dimensions': 'one', 'units': 'meter', 'data': vertical_options["hc"]},
            {'name': 'sc_r', 'long_name': 'S-coordinate at RHO-points', 'datatype': 'f',
             'dimensions': 's_rho', 'units': '-', 'data': self.sc},
            {'name': 'Cs_r', 'long_name': 'S-coordinate stretching curve at RHO-points', 'datatype': 'f',
             'dimensions': 's_rho', 'units': '-', 'data': self.cs}
        ]

    def make(self):
        if not self.var_settings.check_flags():
            return
        # dict of all variables produced
        all_vars = {var['out'] for sublist in [x['variables'] for x in self.sources] for var in sublist}

        # For each group of variables
        for group, group_index in zip(self.sources, count()):
            group_files = ','.join(group['files'])
            weight = group['weight']
            variables = group['variables']
            # set of all variables present in out_file after interpolation
            out_variables = {i["out"] for i in variables}
            # For each file associated with the group of variables
            for in_file, file_index in zip(group['files'], count()):
                # Get the unique filename for each file
                out_file = self.outfile[:-3] + '_' + str(group_index) + '_' + str(file_index) + '.nc'
                # Interpolate everything necessary
                self.interpolator.interpolate(in_file, out_file, weight, variables, group_files)
                # Make any adjustments to variables necessary
                self.var_settings.make_adjustments(out_file, out_variables, all_vars, group_files)
                # Rename time to starting with an underscore if necessary
                if self.time_underscored:
                    self.__rename_time(out_file, group_index)
        # Create a file containing the 1D information
        if has_vertical(self.sources):
            out_file = self.outfile[:-3] + "_" + str(group_index + 1) + "_0.nc"
            self.make_1d(out_file)

    @staticmethod
    def __rename_time(file, group_index):
        with netCDF4.Dataset(file, mode='r+') as my_file:
            my_file.renameDimension('time', group_index + '_time')
            my_file.renameVariable('time', group_index + '_time')

    def make_1d(self, file_name):
        with netCDF4.Dataset(file_name, mode="w")as my_file:  # with automatically opens and closes
            my_file.createDimension('one', 1)
            my_file.createDimension('s_rho', size=len(self.sc))
            for var in [x for x in self.vertical_information if x['name'] not in my_file.variables]:
                file_var = my_file.createVariable(varname=var['name'], datatype=var['datatype'],
                                                  dimensions=(var['dimensions'],))
                file_var[:] = var['data']
                file_var.long_name = var['long_name']
                file_var.units = var['units']

    def add_time_underscores(self):
        self.time_underscored = True
