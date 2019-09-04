from romspy.settings import forcing

target = '/net/kryo/work/munnich/roms/inputs/gridfiles/pactcs30_grd.nc'
out = 'bulk/pactcs30_test_con_frc.nc'
sources = [
    {
        'variables': [
            {'out': 'SST', 'in': 'sst'},
            {'out': 't_air', 'in': 'sat'},
            {'out': 'rho_air', 'in': 'airdens'},
            {'out': 'u_air', 'in': 'w3'},
            {'out': 'humidity', 'in': 'qsea', 'vertical': True}
        ],
        'files': ['/net/kryo/work/updata/ecmwf-reanalysis/era-interim/ERAinterim_flux_landfill.nc'],
        'interpolation_method': 'bil',
    },
    {
        'variables': [
            {'out': 'swrad', 'in': 'shortrad'},
            {'out': 'shflux', 'in': 'netheat'},
            {'out': 'swflux', 'in': 'emp'}
        ],
        'files': ['/net/kryo/work/updata/ecmwf-reanalysis/era-interim/ERAinterim_flux_landfill.nc'],
        'interpolation_method': 'con',
    },
    {
        'variables': [
            {'out': 'sss', 'in': 'salinity'}
        ],
        'files': ['/net/kryo/work/updata/Roms_tools_data/COADS05/sss_landfill.cdf'],
        'interpolation_method': 'bil',
    },
    # {
    #     'variables': [('dust', 'DSTSF')],
    #     'files': ['/net/kryo/work/updata/roms_tools_data/dust/dst79gnx_gx3v5.nc_inputCCSM_fill.nc'],
    #     'interpolation_method': 'bil'
    # },
    {
        'variables': [
            # TODO: Clear up how to do this rotation stuff.
            #  If it's too annoying after all do turning the way you did before!
            {'out': 'tau_y', 'in': 'tauy'},
            {'out': 'tau_x', 'in': 'taux'}
        ],
        'files': ['/work/nicomuen/ERAinterim_tau_landfill.nc'],
        'interpolation_method': 'bil',
    }
]

forcing_object = forcing.Forcing(target, out, sources, verbose=True)
forcing_object.make()

"""
list of dictionaries
            dictionaries have following keys:
             * variables (tuple of variable strings and pairs of strings)
                where pairs of strings if target file has different name for a variable
             * files which is a string of a directory or a list of files
             * interpolation_method which describes the interpolation method


sustr - sms_time = ewss with offset and adjustments
svstr - sms_time = nsss with offset and adjustments
shflux - shf_time = ssr + sshf + slhf + str
pco2_air - pco2_time - opt
swrad - srf_time = ssr
swflux - swf_time = e + tp
precip - swf_time - opt = tp
seaice - ice_time - opt = siconc
SST - sst_time = sst
SSS - sss_time
dQdSST - sss_time = calculation with sst sat rho_atm U qsea
dust - dust_time - opt
iron - iron_time - opt = dust * 62668
"""
