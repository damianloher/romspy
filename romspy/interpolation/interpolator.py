from .shift_grid import adjust_vectors
from .horizontal import calculate_weights, cdo_interpolate
from .vertical import vert_interpolate
from romspy.interpolation.vertical.load_c_libs import bil_weight_extra_len
from romspy.interpolation.vertical import gen_vert_bil, interp_bil
import os
import sys

"""
Author: Nicolas Munnich
License: GNU GPL2+
"""


class Interpolator:
    """
    Interpolates files horizontally and vertically, and can rotate+shift variable pairs if necessary
    Attributes which are not taken from __init__ arguments:
        shift_pairs: ShiftPairCollection object containing variable pairs to shift and rotate if found
    """

    def __init__(self, cdo, target_dir: str, sources: list, target_grid: str, scrip_grid: str, z_levels: tuple,
                 shift_pairs, options: str, in_file: str = '', keep_weights: bool = False,
                 keep_z_clim: bool = False, use_ROMS_grdfile=False, timavg: int = 0, fillmiss=False,
                 lsm_file: str = "", lsm_var: str = "", lsm_limit: float = 0.0, verbose: int = 1):
        """
        Interpolates files horizontally and vertically, and can rotate+shift variable pairs if necessary
        :param cdo: cdo object
        :param sources: list of sources provided
        :param target_grid: ROMS grid file to interpolate onto
        :param scrip_grid: target_grid in SCRIP format
        :param z_levels: 3D array of depths at each point
        :param options: cdo options
        :param in_file: file to interpolate (only needed if sources has no keypdb 'files')
        :param keep_weights: whether to keep weights after program has finished running
        :param keep_z_clim: whether to save the state before vertical interpolation
        :param use_ROMS_grdfile use the ROMS grid file instead of a SCRIP conforming version
        :param timavg number of time records to average over before any interpolations are done
        :param verbose: whether to output runtime information
        """
        self.cdo, self.sources, self.target_grid, self.scrip_grid = cdo, sources, target_grid, scrip_grid
        self.z_levels, self.options, self.keep_weights, self.keep_z_clim = z_levels, options, keep_weights, keep_z_clim
        self.verbose = verbose
        self.in_file = in_file
        self.use_ROMS_grdfile = use_ROMS_grdfile
        self.timavg = timavg
        self.fillmiss = fillmiss
        self.weight_dir = os.path.join(target_dir, "weights")
        os.makedirs(os.path.join(target_dir, "weights"), exist_ok=True)
        self.shift_pairs = shift_pairs
        self.lsm_file = lsm_file
        self.lsm_var = lsm_var
        self.lsm_limit = lsm_limit
        self.outdir = target_dir
        #self.calculate_horizontal_weights()

    def __enter__(self):
        return self

    def __exit__(self, exc_type, exc_val, exc_tb):
        self.clear_weights()

    def interpolate(self, file: str, outfile_name: str, group: dict, group_index: int, variables: list, all_files: str):
        """
        Interpolate a file horizontally, also rotate and shift and vertically interpolate if necessary
        :param file: file to interpolate
        :param outfile_name: output filename
        :param variables: the variables to interpolate from the file
        :param all_files: list of all files containing the variables in variables
        :return:
        """
        if self.verbose > 0:
            print("\033[1;32mInterpolate file horizontally:\033[0m " + file)
        if self.verbose > 1:
            print("Vector pairs saved: " + str(self.shift_pairs.shifts))
        sys.stdout.flush()
        shift_variables = self.shift_pairs.get_shifts(variables)
        shifts = len(shift_variables) > 0
        #vertical_variables = [x['out'] for x in variables if x.get("vertical") is not None]
        vertical_variables = [x['out'] for x in variables if x.get("vertical")]
        vertical = len(vertical_variables) > 0
        if shifts:
            shift_vert_u = self.shift_pairs.get_us(vertical_variables)
            shift_vert_v = self.shift_pairs.get_vs(vertical_variables)
            vertical_variables = [x for x in vertical_variables if x not in shift_vert_u and x not in shift_vert_v]
        if self.verbose > 1:
            print("All Variables: " + str(variables))
            if shifts:
                print("Vectors: " + str(shift_variables))
            if vertical:
                print("Vertical Variables: " + str(vertical_variables))
                if shifts:
                    print(shift_vert_u)
                    print(shift_vert_v)

        # Do time averaging, if necessary:
        if self.timavg > 1:
            import netCDF4
            nc = netCDF4.Dataset(file,'r')
            if 'time' in nc.dimensions:
                nt = len(nc.dimensions['time'])
            elif 'T' in nc.dimensions:
                nt = len(nc.dimensions['T'])
            else:
                msg = 'no time dimension found in: ' + file
                raise ValueError(msg)
            nc.close()
            if nt%self.timavg != 0:
                msg = f'found {nt} time records in: {file}'
                msg += f'averaging over {self.timavg} time records not possible ({self.timavg} does not divide {nt})'
                raise ValueError(msg)
            outfiles = []
            for i in range(int(nt/self.timavg)):
                t1 = i*self.timavg + 1
                t2 = t1 + self.timavg - 1
                input = f'-seltimestep,{t1}/{t2} {file}'
                outfile = self.cdo.timavg(input=input, options=self.options)
                if self.verbose > 0:
                    pdone = f"{int(100*float(t2)/nt):3}% done"
                    if t1 == 1:
                        print(f"   {pdone}", end="")
                        sys.stdout.flush()
                    else:
                        print("\b"*9 + f"{pdone}", end="")
                        sys.stdout.flush()
                outfiles.append(outfile)
            print("")
            # Merge time averages into one file:
            outfile = self.cdo.mergetime(input=(' '.join(outfiles)), options=self.options)
            file = outfile

        if self.keep_z_clim:
            file_path_split = os.path.split(file)
            z_clim_name = os.path.join(file_path_split[0], "z_clim_" + file_path_split[1])

        # Interpolate the file horizontally
        if self.use_ROMS_grdfile:
            grdfile = self.target_grid
        else:
            grdfile = self.scrip_grid
        # Calculate interpolation weights, if necessary:
        self.calculate_horizontal_weights(group, group_index, variables)
        print("weight file used: "+group['weight'])
        sys.stdout.flush()
        outfile = cdo_interpolate(self.cdo, file, group['weight'], grdfile, variables, all_files, self.options,
                                  self.outdir, outfile_name if not (shifts or vertical) else (
                                      z_clim_name if not shifts and self.keep_z_clim else None),
                                  self.lsm_file, self.lsm_var, self.lsm_limit, self.verbose)
        if self.fillmiss:
            # Fill in missing values:
            #outfile = self.cdo.fillmiss2(input=outfile, options=self.options)
            outfile2 = f"{outfile}_2"
            cmd = f"/usr/local/bin/cdo {self.options} -fillmiss2 {outfile} {outfile2}"
            print(cmd)
            sys.stdout.flush()
            os.system(cmd)
            outfile = outfile2

        # Turn and shift variables in shifts
        if shifts:
            if self.verbose > 0:
                print("Shifting and rotating variables: " + str(shift_variables))
            shift_variables_in = []
            for lpair in shift_variables:
                lpair_in = []
                for vdict in variables:
                    if vdict["out"] == lpair[0] or vdict["out"] == lpair[1]:
                        lpair_in.append(vdict["in"])
                shift_variables_in.append(lpair_in)
            outfile = adjust_vectors(self.cdo, outfile, self.target_grid, shift_variables_in, self.options,
                                     self.verbose,
                                     outfile_name if not vertical else (z_clim_name if self.keep_z_clim else None))

        # Interpolate the file vertically
        if vertical:
            if self.verbose > 0:
                print("Interpolating vertically: " + str(vertical_variables))
                sys.stdout.flush()
            if len(vertical_variables) > 0:
                outfile = vert_interpolate(self.cdo, gen_vert_bil, interp_bil, bil_weight_extra_len, outfile,
                                           (outfile_name if not shifts else None), self.weight_dir, vertical_variables,
                                           self.z_levels[0], group, self.options, self.verbose)
            if shifts:
                # noinspection PyTypeChecker
                outfile = vert_interpolate(self.cdo, gen_vert_bil, interp_bil, bil_weight_extra_len, outfile,
                                           None, self.weight_dir, shift_vert_u, self.z_levels[1],
                                           group, self.options, self.verbose)
                outfile = vert_interpolate(self.cdo, gen_vert_bil, interp_bil, bil_weight_extra_len, outfile,
                                           outfile_name, self.weight_dir, shift_vert_v, self.z_levels[2],
                                           group, self.options, self.verbose)
        self.cdo.cleanTempDir()
        return outfile

    def add_shift_pair(self, u: str, v: str):
        """
        Add a pair of variables which should be rotated and shifted
        onto xi_u based grid and eta_v based grid respectively
        :param u: will be rotated and shifted onto xi_u based grid
        :param v: will be rotated and shifted onto eta_v based grid
        :return: None
        """
        self.shift_pairs.add_shift_pair(u, v)

    def calculate_horizontal_weights(self, group, group_index, variables):
        if self.use_ROMS_grdfile:
            grdfile = self.target_grid
        else:
            msg = "use of SCRIP grid file not implemented. Please use the ROMS grid file instead."
            raise ValueError(msg)
        # Check if the source grid is compatible with the weight file:
        import netCDF4
        print("calculate_horizontal_weights: check weight file")
        if "files" in group:
            lfile = group["files"][0]
        else:
            lfile = group["variables"][0]["files"][0]
        #print(f"   lfile = {lfile}")
        nc_grd = netCDF4.Dataset(grdfile,'r')
        nc = netCDF4.Dataset(lfile,'r')
        lvar = variables[0]['in']
        mtd_name = group['interpolation_method']
        weight_name = os.path.join(self.weight_dir, mtd_name + '_weight_g' + str(group_index) + ".nc")
        print(f"   weight_name = {weight_name}")
        sys.stdout.flush()
        if lvar in nc.variables:
            # Check compatibility with source grid (i.e. the data grid):
            vobj = nc.variables[lvar]
            nx = len(nc.dimensions[vobj.dimensions[-1]])
            ny = len(nc.dimensions[vobj.dimensions[-2]])
            if os.path.exists(weight_name):
                nc_w = netCDF4.Dataset(weight_name,'r')
                ns = len(nc_w.dimensions["src_grid_size"])
                nd = len(nc_w.dimensions["dst_grid_size"])
                nc_w.close()
                if ns != nsx*nsy:
                    # Remove weight file, since it was created for a different source grid:
                    print("   remove weight file since it is incompatible with the current source grid:")
                    print(f"      {wname}")
                else:
                    if nd != ndx*ndy:
                        # Remove weight file, since it was created for a different dest grid:
                        print("   remove weight file since it is incompatible with the current dest grid:")
                        os.system(f"rm {wname}")
                        print(f"      {wname}")
                    else:
                        print("   keep weight file since it is compatible with the current source and dest grid:")
                        print(f"      {wname}")
        else:
            print(f"   variable {lvar} not found in lfile -> remove weight file")
            os.system(f"rm {weight_name}")
        nc.close()
        nc_grd.close()
        calculate_weights(self.cdo, self.weight_dir, group, group_index, grdfile, self.options,
                          self.verbose, in_file=self.in_file)

    def clear_weights(self):
        """
        If the weights should be discarded then delete them.
        Otherwise output the information about vertical weights to the console
        :return:
        """
        if not self.keep_weights:
            for group in self.sources:
                try:
                    os.remove(group['weight'])
                except Exception:
                    pass
                try:
                    os.remove(group['vertical_weight'])
                except Exception:
                    pass


class ShiftPairCollection:
    def __init__(self):
        self.shifts = []

    def __contains__(self, item):
        return item in [item for sublist in self.shifts for item in sublist]

    def add_shift_pair(self, x_in: str, y_in: str):
        """
        Add a pair to be shifted. Variables cannot be present in more than one pair, so check for this.
        """
        if x_in not in self and y_in not in self:
            self.shifts.append((x_in, y_in))
        else:
            raise ValueError("ERROR: One of these vectors has already been named in another vector! "
                             "The offending pair: (" + x_in + "," + y_in + ")")

    def get_shifts(self, variables: list):
        """
        Get the pairs which are present in variables.
        """
        out = [x['out'] for x in variables]
        pairs = [(x, y) for x, y in self.shifts if x in out]
        return pairs

    def get_us(self, variables: list):
        us = [x[0] for x in self.shifts]
        return [x for x in us if x in variables]

    def get_vs(self, variables: list):
        vs = [x[1] for x in self.shifts]
        return [x for x in vs if x in variables]
