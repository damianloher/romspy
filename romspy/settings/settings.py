"""
Author: Nicolas Munnich
License: GNU GPL2+
"""


class Settings:
    def __init__(self, adjustments, flags=None):
        self.adjustments = adjustments
        self.flags = {} if flags is None else flags
        self.preprocessor = None

    def make_adjustments(self, file: str, out_variables: set, group_files: str, all_vars: set):  # THIS IS BUGGED
        assert self.preprocessor is not None
        if self.preprocessor.verbose:
            print("Making adjustments to file contents as per settings.")
        for adjustment in self.adjustments:
            if (len(adjustment['out_var_names'] & all_vars) > 0) or (len(adjustment['out_var_names']) == 0):
                contents = adjustment['in_var_names'] & out_variables
                if len(contents) > 0:
                    if len(contents) == len(adjustment['in_var_names']):
                        # If the output isn't already pre-calculated and all the inputs are in the same file
                        adjustment['func'](self.preprocessor, file, group_files, self.flags)
                    else:  # BUGGED
                        print("ERROR: The variables needed to calculate " + str(adjustment['out_var_names']) +
                              " were not all present in the same file! Variables needed: " +
                              str(adjustment['in_var_names']))

    def check_flags(self):
        if self.preprocessor.verbose:
            print("Checking if all adjustments have any necessary flags")
        has_flags = True
        for adjustment in self.adjustments:
            flags = adjustment.get('flags', None)
            if isinstance(flags, list):
                for flag in flags:
                    if flag not in self.flags:
                        print(str(adjustment['func']) + " is missing a necessary flag: " + str(flag))
                        has_flags = False
            elif isinstance(flags, str):
                if flags not in self.flags:
                    print(str(adjustment['func']) + " is missing a necessary flag: " + str(flags))
            else:
                pass
        return has_flags

    def set_flag(self, flag, item):
        self.flags[flag] = item

    def get_flag(self, flag):
        return self.flags[flag]
