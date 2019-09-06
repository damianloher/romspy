from romspy.adjustments import forcing

target = '/net/kryo/work/munnich/roms/inputs/gridfiles/pactcs30_grd.nc'
# scrip = 'ustc90_grd_scripm.nc'
# out = '/net/kryo/work/munnich/Roms/inputs/pacific/pactcs60/pactcs60_test_frc.nc'
out = 'bulk/ERA5_2001_01_on_None.nc'
# sources={'sst':['/net/kryo/work/updata/reynolds_sst/avhrr-only-v2.2002.nc']
sources = [{'variables': [('swrad', 'ssr'), 'ewss', 'nsss', 'sshf', 'slhf', 'str', ('seaice', 'siconc'), ('SST', 'sst'),
                          'e', ('precip', 'tp')],
            'files': ['/net/kryo/work/updata/ecmwf-reanalysis/era5_netcdf/hourly/2001/ERA5_2001_01.nc'],
            'interpolation_method': 'bil'}]

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
