from romspy import PreProcessor, clim_settings

target = '/home/nicomuen/pactcs30_grd.nc'
out = '/home/nicomuen/work/pactcs30_test_clim.nc'
sources = [
    {
        'variables': [
            {'out': 'temp', 'in': 'temp', 'vertical': True},
            {'out': 'salt', 'in': 'salt', 'vertical': True}
        ],
        'files': ['/net/kryo/work/updata/soda_clim/soda_2.1.6/SODA_2.1.6_1979-2008_clm_landfill.nc'],
        'interpolation_method': 'bil',
    }
]

processor = PreProcessor(target, out, sources, clim_settings, keep_weights=True, keep_z_clim=True, verbose=True)
processor.interpolator.add_shift_pair("u", "v")
clim_settings.flags['obc'] = [1, 0, 1, 1]

processor.make()
