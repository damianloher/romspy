import romspy.Global_settings as Global_settings

# How to split the total Chl into SP, Diat, and Diaz:
CHL_SP_RATIO  = 0.9
CHL_DIAT_RATIO = 0.09
CHL_DIAZ_RATIO = 0.01
# Minumum Chla concentration:
CHLAmin = 0.01
# File containing DIC annual mean obtained from CESM (hindcast run 1751-2023):
CESM_DIC_ANN_MEANS = f'{Global_settings.Ocean_3d_input}/CESM_DIC_ann_1751-2023_1x1_landfill.nc'
# First and last year for which CESM_DIC_ANN_MEANS conatins data:
CESM_DIC_ANN_FIRST_YR = 1751
CESM_DIC_ANN_LAST_YR = 2023
