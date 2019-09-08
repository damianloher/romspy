import cdo
import netCDF4
import os
import numpy as np
from itertools import count

from romspy.interpolation.interpolator import ShiftPairCollection
from romspy.verification import test_cdo, has_vertical, verify_sources
from .grid_routines import scrip_grid_from_nc
from .interpolation.vertical import sigma_stretch_sc, sigma_stretch_cs, get_z_levels
from .interpolation import Interpolator

"""
Author: Nicolas Munnich
License: GNU GPL2+
"""


class PreProcessor:
    def __init__(self, target_grid: str, outfile: str, sources: list, **kwargs):
        """
        Contains universal methods

        :param target_grid: grid_routines to be interpolated onto
        :param sources: Information about which variables should be taken from which files using
                which interpolation method onto which type of grid_routines and whether to interpolate them vertically
        :param scrip_grid: Scrip file on rho grid_routines to interpolate with. Will be created if not provided
        :param verbose: whether text should be printed as the program is running
        Can have any of the following optional arguments:
            theta_s - S-coordinate surface control parameter - default 7.0
            theta_b - S-coordinate bottom control parameter - default 0.0
            layers - Number of S-coordinate layers - default 32
            hc - S-coordinate critical depth - default 150
            tcline - S-coordinate surface/bottom layer width - default 150
            sigma_type - default 3
            zeta_source - dict with keys 'name' 'file'. Sets zeta to the first timestep of 'name' in 'file'.
            file_type - output filetype -  default nc4c
            processes - number of processes cdo should use - default 8
            scrip_grid - SCRIP version of target_grid - optional, will be created if not passed
            verbose - whether to print runtime information - default false
            time_underscored - whether 'time' variables should be replaced with '#_time' - default False
            keep_weights - whether to keep calculated weights - default False
            keep_z_clim - whether to keep zclim - default False
        """
        # cdo
        self.cdo = cdo.Cdo()
        test_cdo(self.cdo)
        # Sources
        verify_sources(sources, kwargs.get('verbose', False))
        self.sources, self.target_grid = sources, target_grid

        # Vertical interpolation information
        # Replace anything not passed in with default values
        self.has_vertical = has_vertical(self.sources)
        if self.has_vertical:
            self.theta_s, self.theta_b, self.layers, self.hc, self.tcline, self.sigma_type = (
                kwargs.get('theta_s', 7.0), kwargs.get('theta_b', 0.0),
                kwargs.get("layers", 32), kwargs.get("hc", 150),
                kwargs.get('tcline', 150), kwargs.get("sigma_type", 3)
            )
            self.sc = sigma_stretch_sc(self.layers, True)
            self.cs = sigma_stretch_cs(self.theta_s, self.theta_b, self.sc, self.sigma_type)

            # Get z_levels
            with netCDF4.Dataset(target_grid, mode='r') as my_grid:
                self.h = my_grid.variables['h'][:]
            self.zeta = self.get_zeta(kwargs['zeta_source']) if 'zeta_source' in kwargs else np.zeros_like(self.h)
            self.z_level_rho, self.z_level_u, self.z_level_v = get_z_levels(self.h, self.sc, self.cs, self.hc,
                                                                            self.zeta)

        # CDO options
        self.file_type, self.processes = kwargs.get('file_type', 'nc4c'), kwargs.get('processes', 8)
        # Other Options
        self.verbose, self.time_underscored, self.keep_weights, self.keep_z_clim = (
            kwargs.get('verbose', False), kwargs.get('time_underscored', False),
            kwargs.get('keep_weights', False), kwargs.get('keep_z_clim', False)
        )
        self._adjustments = None
        self.outfile = outfile
        # Interpolator
        self.scrip_grid = kwargs.get('scrip_grid', scrip_grid_from_nc(target_grid))
        self.shift_pairs = ShiftPairCollection()
        if self.verbose:
            print("Finished setup")

    @property
    def options(self):
        return " -b F32 -f " + self.file_type + " -P " + str(self.processes)

    @property
    def adjustments(self):
        return self._adjustments

    @adjustments.setter
    def adjustments(self, value: list):
        if not isinstance(value, list):
            print("ERROR: Your adjustments are in an incorrect format. Adjustments must be a list of dictionaries!")
            return
        for adj in value:
            if not isinstance(adj, dict):
                print("ERROR: Your adjustments are in an incorrect format. Adjustments must be a list of dictionaries!")
                print(adj)
                return
            if not isinstance(adj.get('in_var_names'), set):
                print("ERROR: An adjustment has an incorrect 'in_var_name' key. "
                      "All dicts must have the key 'in_var_names' pointing to a set of input variable strings")
                print(adj)
                return
            if not isinstance(adj.get('out_var_names'), set):
                print("ERROR: An adjustment has an incorrect 'out_var_name' key. "
                      "All dicts must have the key 'out_var_names' pointing to a set of output variable strings")
                print(adj)
                return
            if adj.get('func') is None:
                print("ERROR: All dicts must have the key 'func' pointing to a function!")
                print(adj)
                return
        self._adjustments = value

    def make(self):
        if self.adjustments is None:
            print("Please set your adjustments first!")
            return

        with Interpolator(self.cdo, os.path.split(self.outfile)[0], self.sources, self.target_grid,
                          self.scrip_grid,
                          (self.z_level_rho, self.z_level_u, self.z_level_v) if self.has_vertical else None,
                          self.shift_pairs,
                          self.options, self.keep_weights, self.keep_z_clim, self.verbose) as interp:
            # dict of all variables produced
            all_vars = {var['out'] for sublist in [x['variables'] for x in self.sources] for var in sublist}
            # For each group of variables
            for group, group_index in zip(self.sources, count()):
                group_files = ','.join(group['files'])
                variables = group['variables']
                # set of all variables present in out_file after interpolation
                out_variables = {i["out"] for i in variables}
                # For each file associated with the group of variables
                for in_file, file_index in zip(group['files'], count()):
                    # Get the unique filename for each file
                    out_file = self.outfile[:-3] + '_' + str(group_index) + '_' + str(file_index) + '.nc'
                    # Interpolate everything necessary
                    interp.interpolate(in_file, out_file, group, variables, group_files)
                    # Make any adjustments to variables necessary
                    self.make_adjustments(out_file, out_variables, group_files, all_vars)
                    # Rename time to starting with an underscore if necessary
                    if self.time_underscored:
                        self.__rename_time(out_file, group_index)
                    # add 1D information
                    if has_vertical(self.sources):
                        self.add_1d_attrs(out_file)

    @staticmethod
    def __rename_time(file, group_index):
        with netCDF4.Dataset(file, mode='r+') as my_file:
            my_file.renameDimension('time', group_index + '_time')
            my_file.renameVariable('time', group_index + '_time')

    def add_1d_attrs(self, file_name):
        vertical = [
            {'name': 'theta_s', 'long_name': 'S-coordinate surface control parameter', 'datatype': 'f',
             'dimensions': 'one', 'units': '-', 'data': self.theta_s},
            {'name': 'theta_b', 'long_name': 'S-coordinate bottom control parameter', 'datatype': 'f',
             'dimensions': 'one', 'units': '-', 'data': self.theta_b},
            {'name': 'Tcline', 'long_name': 'S-coordinate surface/bottom layer width', 'datatype': 'f',
             'dimensions': 'one', 'units': 'meter', 'data': self.tcline},
            {'name': 'hc', 'long_name': 'S-coordinate critical depth', 'datatype': 'f',
             'dimensions': 'one', 'units': 'meter', 'data': self.hc},
            {'name': 'sc_r', 'long_name': 'S-coordinate at RHO-points', 'datatype': 'f',
             'dimensions': 's_rho', 'units': '-', 'data': self.sc},
            {'name': 'Cs_r', 'long_name': 'S-coordinate stretching curve at RHO-points', 'datatype': 'f',
             'dimensions': 's_rho', 'units': '-', 'data': self.cs}
        ]
        with netCDF4.Dataset(file_name, mode="r+")as my_file:  # with automatically opens and closes
            for var in vertical:
                my_file.setncattr(var['name'], str(var['long_name'] + " := " + str(var['data'])))

    def add_time_underscores(self):
        self.time_underscored = True

    def make_adjustments(self, file: str, out_variables: set, group_files: str, all_vars: set):
        if self.verbose:
            print("Making adjustments to file contents as per adjustments.")
        for adjustment in self.adjustments:
            # If a variable in out_var_names is not preprovided or it doesn't produce variables
            if len(adjustment['out_var_names'] - all_vars) > 0 or len(adjustment['out_var_names']) == 0:
                # If the file has all the necessary inputs
                if len(adjustment['in_var_names'] & out_variables) == len(adjustment['in_var_names']):
                    # If the output isn't already pre-calculated and all the inputs are in the same file
                    try:
                        if self.verbose:
                            print("Calling " + str(adjustment))
                        adjustment['func'](file, group_files=group_files, options=self.options,
                                           adjustments=self.adjustments, **vars(self))
                    except TypeError as missingflag:
                        print("A flag necessary to run an adjustment was missing, so the adjustment was skipped.")
                        print("Culprit(s): " + str(missingflag))
                        print("File: " + file)
                        print("Please review the documentation on flags and the adjustments you have chosen.")
                        print("If the number of culprit arguments seems extremely long,"
                              " maybe you forgot to add **kwargs to your function?")
                        print("If you still wish to perform the adjustment, then you may call the adjustment manually.")
                elif len(adjustment['in_var_names'] & out_variables) >= 1:
                    print("ERROR: The variables needed to calculate " + str(adjustment['out_var_names']) +
                          " were not all present in the same file! Variables needed: " +
                          str(adjustment['in_var_names']))

    def get_zeta(self, zeta_source):
        temp = self.cdo.remapbil(self.scrip_grid, input="-selname," + zeta_source['name'] + " " + zeta_source['file'],
                                 options=self.options)
        with netCDF4.Dataset(temp) as zetafile:
            zeta = np.array(zetafile.variables[zeta_source['name']][0])
        return zeta

    def mark_as_vectors(self, vector_u_name, vector_v_name):
        self.shift_pairs.add_shift_pair(vector_u_name, vector_v_name)
