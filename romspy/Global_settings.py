use_container = False

if use_container:
    UP_data_dir = '/workspace/datasets'
    era_path = UP_data_dir + '/era5'
    DFS_data_dir = UP_data_dir + '/dfs5.2_factors/'
    Dai_river_runoff_dir = UP_data_dir + '/river_runoff_dai/'
    Ocean_3d_input = UP_data_dir + '/romspy_input'
    ROMS_grd = '/workspace/ROMS'
    ROMS_files = '/workspace/ROMS'
    cdo = "/workspace/.pixi/envs/default/bin/cdo"
    ncks = "/workspace/.pixi/envs/default/bin/ncks"
    ncap2 = "/workspace/.pixi/envs/default/bin/ncap2"
    ncpdq = "/workspace/.pixi/envs/default/bin/ncpdq"
    nccopy = "/workspace/.pixi/envs/default/bin/nccopy"
    ncrename = "/workspace/.pixi/envs/default/bin/ncrename"
else:
    UP_data_dir = '/net/sea/work/datasets'
    era_path = UP_data_dir + '/gridded/atmosphere/2d/reanalysis/era5'
    DFS_data_dir = UP_data_dir + '/gridded/atmosphere/2d/reanalysis/drakkar_forcing_set_5.2/dfs5.2_factors/'
    Dai_river_runoff_dir = UP_data_dir + '/ungridded/surface/ocean/freshwater/river_runoff_dai/'
    Ocean_3d_input = '/net/sea/work/datasets/gridded/ocean/3d/model/romspy_input'
    ROMS_grd = "/net/sea/work/jahaerri/roms/inputs/humpac15_Ncycle/grd"
    ROMS_files = "/net/sea/work/koehne/roms/inputs/humpac15/hindcast_1979_2019/frc/corr_fields_seaice/365days"
    cdo = "/usr/local/bin/cdo"
    ncks = "/usr/local/bin/ncks"
    ncap2 = "/usr/local/bin/ncap2"
    ncpdq = "/usr/local/bin/ncpdq"
    nccopy = "/usr/bin/nccopy"
    ncrename = "/usr/local/bin/ncrename"

cdo_omp_nthreads = 32