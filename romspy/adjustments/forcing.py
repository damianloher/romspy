"""
Author: Nicolas Munnich
License: GNU GPL2+
"""


# TODO: Fix this so that it doesn't use cdo.copy this was only used because critically running out of time

def str_adjustment(file: str, group_files: str, cdo, options, verbose, **kwargs):
    """
    Add attributes to sustr and svstr
    :param preprocessor: PreProcessor object, contains a lot of necessary adjustments
    :param file: file which contains all necessary parameters to calculate dQdSST
    :param group_files: list of all files which contain necessary parameters to calculate dQdSST
    :param flags: any extra info necessary
    :return:
    """
    if verbose:
        print('Adjusting for surface stress')
    temp = cdo.setattribute(
        "sustr@long_name='surface u-momentum stress',svstr@long_name='surface v-momentum stress',"
        "sustr@unit='Newton meter-2',svstr@unit='Newton meter-2',sustr@files='"
        + group_files + "',svstr@files='" + group_files + "'", input=file, options=options)
    cdo.copy(input=temp, output=file, options=options)

    if verbose:
        print('surface stress has been added')


def dust_adjustment(file: str, group_files: str, cdo, options, verbose, **kwargs):
    """
    Calculate iron
    :param preprocessor: PreProcessor object, contains a lot of necessary adjustments
    :param file: file which contains all necessary parameters to calculate dQdSST
    :param group_files: list of all files which contain necessary parameters to calculate dQdSST
    :param flags: any extra info necessary
    :return:
    """
    if verbose:
        print('Adjusting dust to add iron')

    iron_factor = 62668.0
    temp = cdo.setattribute(
        "iron@long_name='iron_flux',iron@units='nmol/cm2/s',iron@files='" + group_files + "'",
        input="-aexpr,'iron=" + str(iron_factor) + " * dust' " + file,
        options=options)
    cdo.copy(input=temp, output=file, options=options)

    if verbose:
        print('Iron has been added')


def swflux_adjustment(file: str, group_files: str, cdo, options, verbose, include_precip=False, **kwargs):
    """
    Calculate surface fresh water flux
    :param preprocessor: PreProcessor object, contains a lot of necessary adjustments
    :param file: file which contains all necessary parameters to calculate dQdSST
    :param group_files: list of all files which contain necessary parameters to calculate dQdSST
    :param flags: any extra info necessary
    :return:
    """
    if verbose:
        print("Adjusting for swflux")
    temp = cdo.setattribute(
        "swflux@long_name='surface fresh water flux',swflux@unit='cm day-1',swflux@files='" + group_files + "'",
        input="-delname,evap" + (
            # 'flags["include_precip"]' instead of 'flags.get("include_precip", False)' makes precip flag mandatory
            ",precip" if not include_precip else "") + " -aexpr,'swflux=precip+evap' " + file,
        options=options)
    cdo.copy(input=temp, output=file, options=options)

    if verbose:
        print("swflux has been added")


def shflux_adjustment(file: str, group_files: str, cdo, options, verbose, **kwargs):
    """
    Compute the surface net heat flux
    :param preprocessor: PreProcessor object, contains a lot of necessary adjustments
    :param file: file which contains all necessary parameters to calculate dQdSST
    :param group_files: list of all files which contain necessary parameters to calculate dQdSST
    :param flags: any extra info necessary
    :return:
    """
    if verbose:
        print('Adjusting for shflux')
    # aexpr to produce the new variable
    # delname removes unnecessary variables
    # setattribute sets the necessary attributes of new variable
    temp = cdo.setattribute(
        "shflux@long_name='surface net heat flux',shflux@unit='Watt meter-2',shflux@files='"
        + group_files + "'",
        input="-delname,lwrad,sensheat,latheat "
              "-aexpr,'shflux=swrad+sshf+slhf+str' " + file,
        options=options)
    cdo.copy(input=temp, output=file, options=options)

    if verbose:
        print('shflux has been added')


def dqdsst_adjustment(file: str, group_files: str, cdo, options, verbose, **kwargs):
    """
        Compute the kinematic surface net heat flux sensitivity to the sea surface temperature: dQdSST
        Q_model ~ Q + dQdSST * (T_model - SST)
        dQdSST = - 4 * eps * stef * T^3 - rho_atm * Cp * CH * U - rho_atm * CE * L * U * 2353 * ln(10 * q_s / T^2)
        SST: sea surface temperature
        t_air: sea surface atmospheric temperature
        rho_air: atmospheric density
        u_air: wind speed
        humidity: sea level specific humidity
        dQdSST: kinematic surface net heat flux sensitivity to the sea surface temperature (Watts meter-2 Celsius-1)
        :param preprocessor: PreProcessor object, contains a lot of necessary adjustments
        :param file: file which contains all necessary parameters to calculate dQdSST
        :param group_files: list of all files which contain necessary parameters to calculate dQdSST
        :param flags: any extra info necessary
    """

    if verbose:
        print("Calculating dQdSST")

    # Specific heat of atmosphere
    Cp = str(1004.8)
    # Sensible heat transfer coefficient (stable condition)
    Ch = str(0.66E-3)
    # Latent heat transfer coefficient (stable condition)
    Ce = str(1.15E-3)
    # Emissivity coefficient
    # eps = str(0.98)
    # Stefan constant
    stef = str(5.6697E-8)
    # SST (Kelvin)
    sstk_formula = "_sst = SST + 273.15;"
    # Latent heat of vaporisation (J.Kg-1)
    l = "(2.008E6 - 2.3E3 * t_air)"
    # Infrared contribution
    q1 = "(-4.0 * " + stef + "* _sst * _sst * _sst)"
    # Sensible heat contribution
    q2 = "(-rho_air * " + Cp + " * " + Ch + "* u_air)"
    # Latent heat contribution
    # 2.30258509299404590109361379290930926799774169921875 = ln(10)
    # 5417.982723814990094979293644428253173828125 = ln(10) * 2353
    dqsdt = "(5417.9827238149900949792936444282531738281250 * humidity / (_sst * _sst))"
    q3 = "(-rho_air * " + Ce + " * " + l + " * u_air * " + dqsdt + " )"

    dqdsst_formula = "'" + sstk_formula + "dQdSST=" + q1 + " + " + q2 + " + " + q3 + "'"

    temp = cdo.setattribute(
        "dQdSST@long_name='surface net heat flux sensitivity to SST',"
        "dQdSST@unit='Watts meter-2 Celsius-1',dQdSST@files='" + group_files + "'",
        input="-delname,sat,rho_atm,U,qsea -aexpr," + dqdsst_formula + " " + file,
        options=options)
    cdo.copy(input=temp, output=file, options=options)

    if verbose:
        print("dQdSST has been added.")


forcing_adjustments = [
    {
        'out_var_names': set(), 'in_var_names': {'sustr', 'svstr'},
        'func': str_adjustment
    },
    {
        'out_var_names': {'shflux'}, 'in_var_names': {'swrad', 'lwrad', 'sensheat', 'latheat'},
        'func': shflux_adjustment
    },
    {
        'out_var_names': {'swflux'}, 'in_var_names': {'evap', 'precip'}, 'func': swflux_adjustment
    },
    {
        'out_var_names': {'iron'}, 'in_var_names': {'dust'}, 'func': dust_adjustment
    },
    {
        'out_var_names': {'dQdSST'}, 'in_var_names': {'SST', 't_air', 'rho_air', 'u_air', 'humidity'},
        'func': dqdsst_adjustment
    }
]

#
# dict={
#     'out_var_name':['a'],
#     'in_var_name':['b','c'],
#     'func': b_c_to_a
# }
#
# forcing_adjustments.append(dict)
