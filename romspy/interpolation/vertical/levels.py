import numpy as np
from math import sinh, cosh, tanh, exp
from romspy.interpolation.shift_grid import shift

"""
Author: Nicolas Munnich
License: GNU GPL2+
"""


def sigma_stretch_sc(num_levels: int, rho_grid: bool = True, verbose=False) -> np.ndarray:
    """
    Compute S-coordinate level locations
    :param num_levels: depth layers
    :param rho_grid: True if vertical grid type is rho point, False if vertical grid type is w point
    :param verbose: If true prints a message when this function is called
    :return: S-coordinate level locations
    """
    if verbose:
        print("Calculating S-coord level locations")
    if rho_grid:
        return (np.arange(1, num_levels + 1) - num_levels - 0.5) / num_levels
    else:
        return (np.arange(0, num_levels + 1) - num_levels) / num_levels


def sigma_stretch_cs(theta_s: float, theta_b: float, sc: np.ndarray = None, sigma_type: int = 3, verbose=False):
    """
    Compute S-coordinate stretching factor Cs

    :param theta_s: stretching parameter surface
    :param theta_b: stretching parameter bottom
    :param sc: level locations
    :param sigma_type: sigma coordinate system type.
    :param verbose: If true prints a message when this function is called
    :return: S-coordinate stretching factor
    """

    if verbose:
        print("Calculating S-coordinate stretching factor")

    if sigma_type <= 2:
        cff1 = 1 / sinh(theta_s)
        cff2 = 0.5 / tanh(0.5 * theta_s)
        return (1 - theta_b) * cff1 * np.sinh(sc * theta_s) + theta_b * (cff2 * np.tanh(theta_s * (sc + 0.5)) - 0.5)

    elif sigma_type == 3:
        if theta_s > 0:
            csrf = (1 - np.cosh(theta_s * sc)) / (cosh(theta_s) - 1)
        else:
            csrf = -np.power(sc, 2)

        if theta_b > 0:
            return np.expm1(theta_b * csrf) / (1 - exp(-theta_b))
        else:
            return csrf
    else:
        raise ValueError("Unknown Sigma coord type.")


def z_levels(h: np.ndarray, sc: np.ndarray, cs: np.ndarray, hc: int, zeta: np.ndarray = None, sigma_type=3,
             verbose=False) -> np.ndarray:
    """
    Compute the depth of each point in a 3D grid
    :param h: 2D array of seafloor height at each point
    :param sc: s-coord level locations
    :param cs: s-coord stretching factor
    :param hc: S-coordinate critical depth
    :param zeta: 2D array of sea surface height
    :param sigma_type: s-coord type
    :param verbose:
    :return: depth of points as 3D array shaped sc * h.shape[0] * h.shape[1]
    """
    if zeta is None:
        zeta = np.zeros_like(h)
    if sigma_type == 1:
        if verbose:
            print("Using s-coord type 1 (old)")
            print("h shape is: " + str(h.shape))

        hinv = 1 / h
        cff = hc * (sc - cs)
        z0 = cff[:, None, None] + cs[:, None, None] * h
        z = z0 + zeta * (1 + z0 * hinv)
    elif sigma_type == 2:
        if verbose:
            print("Using s-coord type 2 (ETH 1.0 legacy)")
            print("h shape is: " + str(h.shape))

        hinv = 1 / (h + hc)
        cff = hc * sc
        z = zeta + (zeta + h) * (cff[:, None, None] + cs[:, None, None] * h) * hinv
    elif sigma_type == 3:
        if verbose:
            print("Using s-coord type 3 (new)")
            print("h shape is: " + str(h.shape))
        hinv = 1 / (h + hc)
        cff = hc * sc
        z = zeta + (zeta + h) * (cff[:, None, None] + cs[:, None, None] * h) * hinv
    else:
        raise ValueError("Unknown Sigma coord type.")

    return z


def get_z_levels(h: np.ndarray, sc: np.ndarray, cs: np.ndarray, hc: int, zeta: np.ndarray = None, sigma_type=3,
                 verbose=False) -> tuple:
    """
    Compute the depth of each point in a 3D grid
    :param h: 2D array of seafloor height at each point
    :param sc: s-coord level locations
    :param cs: s-coord stretching factor
    :param hc: S-coordinate critical depth
    :param zeta: 2D array of sea surface height
    :param sigma_type: s-coord type
    :param verbose:
    :return: depth of points as 3D array shaped sc * h.shape[0] * h.shape[1]
    """
    my_z_levels = (z_levels(h, sc, cs, hc, zeta, sigma_type, verbose),
                   z_levels(shift(h, 1), sc, cs, hc, shift(zeta, 1), sigma_type, verbose),
                   z_levels(shift(h, 2), sc, cs, hc, shift(zeta, 2), sigma_type, verbose))

    return my_z_levels
