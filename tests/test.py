import numpy as np
from math import sinh, cosh, tanh, exp
import sys


class Test:

    def __init__(self):
        self.sc = None
        self.cs = None
        self.hc = None
        self._1d_and_scalar_vars = None
        self.add_vertical_options()
        a = np.tile(np.array([5, 5, 5, 5, 5]), (1, 10))
        zlevs = self.zlevs(a)
        print(zlevs)

    def add_vertical_options(self, ):
        self.sc = self.sigma_stretch_sc(32, True)
        self.cs = self.sigma_stretch_cs(7,
                                        0, self.sc, 3)
        self.hc = 150
        self._1d_and_scalar_vars = [
            {'name': 'theta_s', 'long_name': 'S-coordinate surface control parameter', 'datatype': 'f',
             'dimensions': 'one', 'units': '-', 'data': 7},
            {'name': 'theta_b', 'long_name': 'S-coordinate bottom control parameter', 'datatype': 'f',
             'dimensions': 'one', 'units': '-', 'data': 0},
            {'name': 'Tcline', 'long_name': 'S-coordinate surface/bottom layer width', 'datatype': 'f',
             'dimensions': 'one', 'units': 'meter', 'data': 150},
            {'name': 'hc', 'long_name': 'S-coordinate critical depth', 'datatype': 'f',
             'dimensions': 'one', 'units': 'meter', 'data': 150},
            {'name': 'sc_r', 'long_name': 'S-coordinate at RHO-points', 'datatype': 'f',
             'dimensions': 's_rho', 'units': '-', 'data': self.sc},
            {'name': 'Cs_r', 'long_name': 'S-coordinate stretching curve at RHO-points', 'datatype': 'f',
             'dimensions': 's_rho', 'units': '-', 'data': self.cs}
        ]

    @staticmethod
    def sigma_stretch_sc(sc_num: int, rho_type: bool = True) -> np.ndarray:
        """
        Compute S-coordinate level locations
        :param sc_num: depth layers
        :param rho_type: True if grid_routines type is rho point, False if grid_routines type is w point
        :return: S-coordinate level locations
        """
        if rho_type:
            return (np.arange(1, sc_num + 1) - sc_num - 0.5) / sc_num
        else:
            return (np.arange(0, sc_num + 1) - sc_num) / sc_num

    @staticmethod
    def sigma_stretch_cs(theta_s: float, theta_b: float, sc: np.ndarray, sigma_type: int = 3):
        """
        Compute S-coordinate stretching factor Cs

        :param theta_s: stretching parameter surface
        :param theta_b: stretching parameter bottom
        :param sc: level locations
        :param sigma_type: sigma coordinate system type.
        :return: S-coordinate stretching factor
        """

        cff1 = 1 / sinh(theta_s)
        cff2 = 0.5 / tanh(0.5 * theta_s)

        if sigma_type <= 2:
            return (1 - theta_b) * cff1 * np.sinh(theta_s * sc) + theta_b * (cff2 * np.tanh(theta_s * (sc + 0.5)) - 0.5)
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

    def zlevs(self, h: np.ndarray, zeta: np.ndarray = None, sigma_type=3, verbose=False) -> np.ndarray:
        """
        Compute the depth of points
        :param h:2D array of seafloor height at each point
        :param zeta: 2D array of sea surface height
        :param sigma_type:
        :param verbose:
        :return: depth of points as 3D array shaped sc * h.shape[0] * h.shape[1]
        """
        if self.sc is None or self.cs is None or self.hc is None:
            sys.exit("ERROR: Trying to vertically interpolate without sc, cs, hc. What have you done?!")
        if zeta is None:
            zeta = np.zeros_like(h)
        if sigma_type == 1:
            if verbose:
                print("Using s-coord type 1 (old)")
                print("h shape is: " + str(h.shape))

            hinv = 1 / h
            cff = self.hc * (self.sc - self.cs)
            z0 = cff[:, None, None] + self.cs[:, None, None] * h
            z = z0 + zeta * (1 + z0 * hinv)
        elif sigma_type == 2:
            if verbose:
                print("Using s-coord type 2 (ETH 1.0 legacy)")
                print("h shape is: " + str(h.shape))

            hinv = 1 / (h + self.hc)
            cff = self.hc * self.sc
            z = zeta + (zeta + h) * (cff[:, None, None] + self.cs[:, None, None] * h) * hinv
        elif sigma_type == 3:
            if verbose:
                print("Using s-coord type 3 (new)")
                print("h shape is: " + str(h.shape))
            hinv = 1 / (h + self.hc)
            cff = self.hc * self.sc
            z = zeta + (zeta + h) * (cff[:, None, None] + self.cs[:, None, None] * h) * hinv
        else:
            raise ValueError("Unknown Sigma coord type.")

        return z


Test()
