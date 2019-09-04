from romspy.settings import clim

target = '/net/kryo/work/munnich/roms/inputs/gridfiles/pactcs30_grd.nc'
out = 'bulk/pactcs30_test_clim.nc'
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

clim_obj = clim.Climatology(target, out, sources, verbose=True)
clim_obj.make()
