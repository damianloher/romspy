from .shift_grid import shift_rotate_grid
from .horizontal import calculate_weights, cdo_interpolate
from .vertical import vert_interpolate
from romspy.interpolation.vertical.load_c_libs import bil_weight_extra_len
from romspy.interpolation.vertical import gen_vert_bil, interp_bil
import os

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

    def __init__(self, cdo, target_dir: str, sources: list, target_grid: str, scrip_grid: str, z_levels: dict,
                 options: str, vertical_weights: dict = None, keep_weights: bool = False, keep_z_clim: bool = False,
                 verbose: bool = False):
        """
        Interpolates files horizontally and vertically, and can rotate+shift variable pairs if necessary
        :param cdo: cdo object
        :param sources: list of sources provided
        :param target_grid: grid to interpolate onto
        :param scrip_grid: target_grid in SCRIP format
        :param z_levels: 3D array of depths at each point
        :param options: cdo options
        :param vertical_weights: dict which is used to keep track of created vertical weights, to reapply and delete
        :param keep_weights: whether to keep weights after program has finished running
        :param keep_z_clim: whether to save the state before vertical interpolation
        :param verbose: whether to output runtime information
        """
        self.cdo = cdo
        self.sources = sources
        self.target_grid = target_grid
        self.scrip_grid = scrip_grid
        self.z_levels = z_levels
        self.options = options
        self.keep_weights = keep_weights
        self.keep_z_clim = keep_z_clim
        self.verbose = verbose
        self.weight_dir = os.path.join(target_dir, "weights")
        if not os.path.exists(os.path.join(target_dir, "weights")):
            os.mkdir(os.path.join(target_dir, "weights"))
        self.vertical_weights = {} if vertical_weights is None else vertical_weights
        self.shift_pairs = ShiftPairCollection()
        self.calculate_horizontal_weights()

    """"with netCDF4.Dataset(target_grid) as grd:
        h_vals = shift(grd.variables['h'][:], tuple(dim_lens[-2:]))
    if ssh:
        ssh_vals = shift(my_file.variables['ssh'][:], tuple(dim_lens[-2:]))
    z_levels = get_z_levels(h_vals, sc, cs, hc, ssh_vals, sigma_type=3, verbose=verbose)
    """

    def __enter__(self):
        return self

    def __exit__(self, exc_type, exc_val, exc_tb):
        self.clear_weights()

    def interpolate(self, file: str, outfile_name: str, weight: str, variables: list, all_files: str):
        """
        Interpolate a file horizontally, also rotate and shift and vertically interpolate if necessary
        :param file: file to interpolate
        :param outfile_name: output filename
        :param weight: the horizontal weight belonging to the file
        :param variables: the variables to interpolate from the file
        :param all_files: list of all files containing the variables in variables
        :return:
        """
        shift_variables = self.shift_pairs.get_shift_pairs(variables)
        shifts = len(shift_variables) > 0
        vertical_variables = [x['out'] for x in variables if x.get("vertical") is not None]
        vertical = len(vertical_variables) > 0
        shift_vert_u = [x for x in vertical_variables if x in self.shift_pairs.u_list]
        shift_vert_v = [x for x in vertical_variables if x in self.shift_pairs.v_list]
        vertical_variables = [x for x in vertical_variables if (x not in shift_vert_u) and (x not in shift_vert_v)]

        if self.keep_z_clim:
            file_path_split = os.path.split(file)
            z_clim_name = os.path.join(file_path_split[0], "z_clim_" + file_path_split[1])

        # Interpolate the file horizontally
        outfile = cdo_interpolate(self.cdo, file, weight, self.scrip_grid, variables, all_files, self.options,
                                  outfile_name if not (shifts or vertical) else (
                                      z_clim_name if not shifts and self.keep_z_clim else None),
                                  self.verbose)

        # Turn and shift variables in shifts
        if shifts:
            outfile = shift_rotate_grid(self.cdo, outfile, self.target_grid, shift_variables, self.options, self.verbose,
                                        outfile_name if not vertical else (z_clim_name if self.keep_z_clim else None))

        # Interpolate the file vertically
        if vertical:
            if len(vertical_variables) > 0:
                outfile = vert_interpolate(self.cdo, gen_vert_bil, interp_bil, bil_weight_extra_len, outfile,
                                           outfile_name, self.weight_dir, vertical_variables, self.z_levels['rho-rho'],
                                           self.vertical_weights, self.options, self.verbose)
            if len(shift_vert_u) > 0:
                outfile = vert_interpolate(self.cdo, gen_vert_bil, interp_bil, bil_weight_extra_len, outfile,
                                           outfile_name, self.weight_dir, shift_vert_u, self.z_levels['rho-u'],
                                           self.vertical_weights, self.options, self.verbose)
            if len(shift_vert_v) > 0:
                outfile = vert_interpolate(self.cdo, gen_vert_bil, interp_bil, bil_weight_extra_len, outfile,
                                           outfile_name, self.weight_dir, shift_vert_v, self.z_levels['v-rho'],
                                           self.vertical_weights, self.options, self.verbose)
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

    def calculate_horizontal_weights(self):
        calculate_weights(self.cdo, self.weight_dir, self.sources, self.scrip_grid, self.options, self.verbose)

    def clear_weights(self):
        """
        If the weights should be discarded then delete them.
        Otherwise output the information about vertical weights to the console
        :return:
        """
        if not self.keep_weights:
            for file, ssh in self.vertical_weights.values():
                os.remove(file)
            for group in self.sources:
                os.remove(group['weight'])
        else:
            print("The vertical weights along with usage information are:")
            print(self.vertical_weights)
            print("Please copy and save if you wish to reuse the weights.")
            print("This should be passed in as the 'vertical_weights' argument.")


class ShiftPairCollection:
    def __init__(self):
        self.u_list = []
        self.v_list = []

    def __contains__(self, item):
        return item in self.u_list or item in self.v_list

    def add_shift_pair(self, u: str, v: str):
        """
        Add a pair to be shifted. Variables cannot be present in more than one pair, so check for this.
        """
        if u not in self and v not in self:
            self.u_list.append(u)
            self.v_list.append(v)
        else:
            raise ValueError("ERROR: Variable flagged for shifting has already been flagged elsewhere! "
                             "The offending pair: (" + u + "," + v + ")")

    def get_shift_pairs(self, variables: list):
        """
        Get the pairs which are present in variables. If one variable is present but not the other, raises a ValueError.
        """
        pairs = [(u, v) for u, v in zip(self.u_list, self.v_list) if u in variables and v in variables]
        for var in variables:
            if var in self.u_list or var in self.v_list:
                zip_pairs = list(zip(pairs))
                if var not in zip_pairs[0] and var not in zip_pairs[1]:
                    raise ValueError(
                        "ERROR: Variable flagged for shifting does not have its pair present! Variable: " + var)
        return pairs
