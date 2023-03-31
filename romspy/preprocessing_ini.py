import cdo
import netCDF4
import os
import numpy as np
from itertools import count

from romspy.interpolation.interpolator import ShiftPairCollection
from romspy.verification import test_cdo, has_vertical, verify_sources
#from .grid_routines import scrip_grid_from_nc
from .interpolation.vertical import sigma_stretch_sc, sigma_stretch_cs, get_z_levels
#from .interpolation import Interpolator

"""
Author: Damian Loher, based on code from Nicolas Munnich
License: GNU GPL2+
"""


class PreProcessorIni:
    def __init__(self, outfile: str, clim_files: list, **kwargs):
        """
        Contains universal methods for generating initial conditions

        :param clim_files: list of paths to clim files to be used for the initial conditions
        :param verbose: whether text should be printed as the program is running
        Can have any of the following optional arguments:
            theta_s - S-coordinate surface control parameter - default 7.0
            theta_b - S-coordinate bottom control parameter - default 0.0
            layers - Number of S-coordinate layers - default 32
            hc - S-coordinate critical depth - default 150
            tcline - S-coordinate surface/bottom layer width - default 150
            sigma_type - default 3
            file_type - output filetype -  default nc4c
            processes - number of processes cdo should use - default 8
            verbose - whether to print runtime information - default false
        """
        self.verbose = kwargs.get('verbose', False)
        if self.verbose:
            print()
            print('\033[1;32m==================\033[0m')
            print('\033[1;32mInitial conditions\033[0m')
            print('\033[1;32m==================\033[0m')
            print()
        # cdo
        self.cdo = cdo.Cdo(debug=kwargs.get('verbose', False))
        test_cdo(self.cdo)
        self.clim_files = []
        for f in clim_files:
            if os.path.exists(f):
                self.clim_files.append(f)
        if len(self.clim_files) == 0:
            raise ValueError('no clim files found')

        # Information about the sigma coordinates:
        # Replace anything not passed in with default values
        self.theta_s, self.theta_b, self.layers, self.hc, self.tcline, self.sigma_type = (
            kwargs.get('theta_s', 7.0), kwargs.get('theta_b', 0.0),
            kwargs.get("layers", 32), kwargs.get("hc", 150),
            kwargs.get('tcline', 150), kwargs.get("sigma_type", 3)
        )
        self.sc = sigma_stretch_sc(self.layers, True)
        self.cs = sigma_stretch_cs(self.theta_s, self.theta_b, self.sc, self.sigma_type)

        # CDO options
        self.file_type, self.processes = kwargs.get('file_type', 'nc4c'), kwargs.get('processes', 8)
        # Other Options
        self._adjustments = None
        self.outfile = outfile
        if self.verbose:
            print("Finished PreProcessorIni setup")

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

        # Loop over clim files and generate corresponding initial conditions:
        fcnt = 1
        for fclm in self.clim_files:
            if not os.path.exists(fclm):
                print("File not found and skipped: "+fclm)
            # Get number of time records in fclm:
            nc = netCDF4.Dataset(fclm,'r')
            nt = len(nc.dimensions['time'])
            nc.close()
            # Determine output file name:
            fdir, fname = os.path.split(fclm)
            if 'clm' in fname:
                outfile = fdir + '/' + fname.replace('clm','ini')
            elif 'clim' in fname:
                outfile = fdir + '/' + fname.replace('clim','ini')
            else:
                tmp = fname.split('_')
                if len(tmp)>=3:
                    outfile = fdir + '/' + tmp[:-2] + '_ini_' + tmp[-2] + '_' + tmp[-1]
                else:
                    outfile = fdir + '/' + 'ini_{:02}.nc'.format(fcnt)
            fcnt += 1
            # Execute cdo command for averaging over the first and last time record of the clim file:
            print("Generate ini file: "+outfile)
            self.cdo.timavg(input=(' -seltimestep,1,{} {}'.format(nt,fclm)), options=self.options, output=outfile)
            self.add_1d_attrs(outfile)
            # Set time variable to 0:
            nc = netCDF4.Dataset(outfile,'a')
            vobj = nc.variables['time']
            vobj[:] = 0.0
            nc.close()

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
        with netCDF4.Dataset(file_name, mode="r+") as my_file:  # with automatically opens and closes
            for var in vertical:
                #my_file.setncattr(var['name'], str(var['long_name'] + " := " + str(var['data'])))
                my_file.setncattr(var['name'], var['data'])

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
