import signal
import sys
import os


def verify_sources(sources, verbose: bool):
    """
    :param sources: list(
                            dict{
                                'variables' : list(
                                                    {
                                                        'in' : var,
                                                        'out' : var
                                                    }
                                                ),
                                'files' : list(file),
                                'interpolation_method' : ('bil' | 'bic' | 'nn' | 'dis' | 'con' | 'con2' | 'laf'),
                                'vertical' : (True | False),
                                'grid_type' : ('rho' | 'u' | 'v')
                                }
                            )
    :param verbose: print runtime info
    """
    if verbose:
        print("Verifying that the sources were input correctly")
    for a_dict in sources:
        try:
            for var_dict in a_dict['variables']:
                try:
                    for item in ['in', 'out']:
                        assert item in var_dict.keys()
                except AssertionError:
                    sys.exit("ERROR: All variables must have 'in' and 'out' values! Offending dict: " + str(var_dict))

            files = a_dict['files']
            try:
                if isinstance(files, list):
                    for i in files:
                        assert os.path.isfile(i)
                elif os.path.isdir(files):
                    a_dict['files'] = get_files_in_dir(files)
                else:
                    assert os.path.isfile(files)
                    a_dict['files'] = [files]

            except AssertionError:
                sys.exit("ERROR: Not all filepaths are valid! Offending file list: " + str(files))

            try:
                assert a_dict['interpolation_method'] is not None
            except KeyError:
                sys.exit("Missing interpolation method!")
        except KeyError:
            sys.exit("ERROR: A dictionary in sources is missing a mandatory key."
                     "\nKeys: 'variables', 'files', 'interpolation_method'")


def get_files_in_dir(root: str) -> list:
    result = []
    for path, subdirs, files in os.walk(root):
        for name in files:
            result.append(os.path.join(path, name))
    return result


def test_cdo(my_cdo):
    """
    Calls cdo help to test if cdo is functioning properly.
    Aborts and informs user if this call takes more than 5 seconds.
    :param my_cdo: cdo object
    :return:
    """

    def cdo_handler(signum, frame):
        raise IOError(
            "CDO likely does not exist! Either that or it is running far too slowly for use. "
            "Try 'module load cdo' or switching to a different computer. Make sure CDO is installed!")

    signal.signal(signal.SIGALRM, cdo_handler)
    signal.alarm(5)
    try:
        my_cdo.sinfov
    except IOError as err:
        print(err)
        sys.exit()
    finally:
        signal.signal(signal.SIGALRM, signal.SIG_IGN)


def has_vertical(sources: list) -> bool:
    """
    Checks if the list of sources has a vertical interpolation anywhere
    :param sources:
    :return:
    """
    for group in sources:
        for variable in group['variables']:
            if variable.get('vertical', False):
                return True
    return False
