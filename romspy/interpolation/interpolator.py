from .shift_grid import adjust_vectors
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

    def __init__(self, cdo, target_dir: str, sources: list, target_grid: str, scrip_grid: str, z_levels: tuple,
                 shift_pairs, options: str, keep_weights: bool = False, keep_z_clim: bool = False,
                 verbose: bool = False):
        """
        Interpolates files horizontally and vertically, and can rotate+shift variable pairs if necessary
        :param cdo: cdo object
        :param sources: list of sources provided
        :param target_grid: grid to interpolate onto
        :param scrip_grid: target_grid in SCRIP format
        :param z_levels: 3D array of depths at each point
        :param options: cdo options
        :param keep_weights: whether to keep weights after program has finished running
        :param keep_z_clim: whether to save the state before vertical interpolation
        :param verbose: whether to output runtime information
        """
        self.cdo, self.sources, self.target_grid, self.scrip_grid = cdo, sources, target_grid, scrip_grid
        self.z_levels, self.options, self.keep_weights, self.keep_z_clim = z_levels, options, keep_weights, keep_z_clim
        self.verbose = verbose
        self.weight_dir = os.path.join(target_dir, "weights")
        if not os.path.exists(os.path.join(target_dir, "weights")):
            os.mkdir(os.path.join(target_dir, "weights"))
        self.shift_pairs = shift_pairs
        self.calculate_horizontal_weights()

    def __enter__(self):
        return self

    def __exit__(self, exc_type, exc_val, exc_tb):
        self.clear_weights()

    def interpolate(self, file: str, outfile_name: str, group: dict, variables: list, all_files: str):
        """
        Interpolate a file horizontally, also rotate and shift and vertically interpolate if necessary
        :param file: file to interpolate
        :param outfile_name: output filename
        :param variables: the variables to interpolate from the file
        :param all_files: list of all files containing the variables in variables
        :return:
        """
        self.shift_pairs: ShiftPairCollection
        shift_variables = self.shift_pairs.get_shifts(variables)
        shifts = len(shift_variables) > 0
        vertical_variables = [x['out'] for x in variables if x.get("vertical") is not None]
        vertical = len(vertical_variables) > 0
        shift_vert_u = self.shift_pairs.get_us(vertical_variables)
        shift_vert_v = self.shift_pairs.get_vs(vertical_variables)
        vertical_variables = list(set(vertical_variables) - (set(shift_vert_u) | set(shift_vert_v)))

        if self.keep_z_clim:
            file_path_split = os.path.split(file)
            z_clim_name = os.path.join(file_path_split[0], "z_clim_" + file_path_split[1])

        if self.verbose:
            print("Interpolating horizontally")
        # Interpolate the file horizontally
        outfile = cdo_interpolate(self.cdo, file, group['weight'], self.scrip_grid, variables, all_files, self.options,
                                  outfile_name if not (shifts or vertical) else (
                                      z_clim_name if not shifts and self.keep_z_clim else None),
                                  self.verbose)

        # Turn and shift variables in shifts
        if shifts:
            if self.verbose:
                print("Shifting and rotating variables: " + str(shift_variables))
            outfile = adjust_vectors(self.cdo, outfile, self.target_grid, shift_variables, self.options,
                                     self.verbose,
                                     outfile_name if not vertical else (z_clim_name if self.keep_z_clim else None))

        # Interpolate the file vertically
        if vertical:
            if self.verbose:
                print("Interpolating vertically: " + str(vertical_variables))
            if len(vertical_variables) > 0:
                outfile = vert_interpolate(self.cdo, gen_vert_bil, interp_bil, bil_weight_extra_len, outfile,
                                           outfile_name, self.weight_dir, vertical_variables, self.z_levels[0],
                                           group, self.options, self.verbose)
            if len(shift_vert_u) > 0:
                outfile = vert_interpolate(self.cdo, gen_vert_bil, interp_bil, bil_weight_extra_len, outfile,
                                           outfile_name, self.weight_dir, shift_vert_u, self.z_levels[1],
                                           group, self.options, self.verbose)
            if len(shift_vert_v) > 0:
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

    def calculate_horizontal_weights(self):
        calculate_weights(self.cdo, self.weight_dir, self.sources, self.scrip_grid, self.options, self.verbose)

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
        self.shifts = {}

    def __contains__(self, item):
        return item in [item for sublist in self.shifts.keys() for item in sublist]

    def add_shift_pair(self, x_in: str, y_in: str, u_out: str, v_out: str):
        """
        Add a pair to be shifted. Variables cannot be present in more than one pair, so check for this.
        """
        if x_in not in self and y_in not in self:
            self.shifts[(x_in, y_in)] = (u_out, v_out)
        else:
            raise ValueError("ERROR: One of these vectors has already been named in another vector! "
                             "The offending pair: (" + x_in + "," + y_in + ")")

    def get_shifts(self, variables: list):
        """
        Get the pairs which are present in variables.
        """
        pairs = [((x, y), self.shifts.get((x, y))) for x, y in zip(variables, variables) if
                 self.shifts.get((x, y)) is not None]
        return pairs

    def get_us(self, variables: list):
        return [x[1][0] for x in self.shifts.items() if x[0][0] in variables]

    def get_vs(self, variables: list):
        return [x[1][1] for x in self.shifts.items() if x[0][1] in variables]
