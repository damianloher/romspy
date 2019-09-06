from romspy import PreProcessor, forcing_adjustments

target = 'pactcs30_grd.nc'
out = 'test_pactcs30_frc.nc'
sources = [
    {
        'variables': [
            {'out': 'SST', 'in': 'sst'},
            {'out': 't_air', 'in': 'sat'},
            {'out': 'rho_air', 'in': 'airdens'},
            {'out': 'u_air', 'in': 'w3'},
            {'out': 'humidity', 'in': 'qsea'},
            {'out': 'swrad', 'in': 'shortrad'},
            {'out': 'shflux', 'in': 'netheat'},
            {'out': 'swflux', 'in': 'emp'}
        ],
        'files': ['/net/kryo/work/updata/ecmwf-reanalysis/era-interim/ERAinterim_flux_landfill.nc'],
        'interpolation_method': 'bil',
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
    # {
    #     'variables': [
    #         {'out': 'sustr', 'in': 'tauy'},
    #         {'out': 'svstr', 'in': 'taux'}
    #     ],
    #     'files': ['/work/nicomuen/ERAinterim_tau_landfill.nc'], # This file got deleted, find a new one
    #     'interpolation_method': 'bil',
    # }
]

processor = PreProcessor(target, out, sources, forcing_settings, keep_weights=True, keep_z_clim=True, verbose=True)
processor.interpolator.add_shift_pair("sustr", "svstr")
forcing_settings.flags['include_precip'] = False

processor.make()
