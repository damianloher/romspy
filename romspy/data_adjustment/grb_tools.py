import sys
import os

"""
Author: Nicolas Munnich
License: GNU GPL2+
"""
# Dict of .grb era interim -> classic variable name and attributes
para_table_era_interim = [
    {'name': 'var0', 'out_name': 'var0', 'long_name': 'undefined', 'units': ''},
    {'name': 'var1', 'out_name': 'None', 'long_name': 'Stream Function ', 'units': 'm**2 s**-1'},
    {'name': 'var2', 'out_name': 'None', 'long_name': 'Velocity Potential', 'units': ' m**2 s**-1'},
    {'name': 'var3', 'out_name': 'None', 'long_name': 'Potential Temperature ', 'units': 'K'},
    {'name': 'var4', 'out_name': 'None', 'long_name': 'Equivalent Potential Temperature', 'units': ' K'},
    {'name': 'var5', 'out_name': 'None', 'long_name': 'Saturated Equivalent Potential Temperature', 'units': ' K'},
    {'name': 'var6', 'out_name': 'None', 'long_name': 'RESERVED FOR METVIEW. ', 'units': '-'},
    {'name': 'var7', 'out_name': 'None', 'long_name': 'RESERVED FOR METVIEW.', 'units': ' -'},
    {'name': 'var8', 'out_name': 'None', 'long_name': 'RESERVED FOR METVIEW. ', 'units': '-'},
    {'name': 'var9', 'out_name': 'None', 'long_name': 'RESERVED FOR METVIEW. ', 'units': '-'},
    {'name': 'var10', 'out_name': 'None', 'long_name': 'RESERVED FOR METVIEW. ', 'units': '-'},
    {'name': 'var11', 'out_name': 'None', 'long_name': 'u-component of Divergent Wind ', 'units': 'm s**-1'},
    {'name': 'var12', 'out_name': 'None', 'long_name': 'v-component of Divergent Wind ', 'units': 'm s**-1'},
    {'name': 'var13', 'out_name': 'None', 'long_name': 'u-component of Rotational Wind ', 'units': 'm s**-1'},
    {'name': 'var14', 'out_name': 'None', 'long_name': 'v-component of Rotational Wind ', 'units': 'm s**-1'},
    {'name': 'var15', 'out_name': 'None', 'long_name': 'RESERVED FOR METVIEW. ', 'units': '-'},
    {'name': 'var16', 'out_name': 'None', 'long_name': 'RESERVED FOR METVIEW. ', 'units': '-'},
    {'name': 'var17', 'out_name': 'None', 'long_name': 'RESERVED FOR METVIEW. ', 'units': '-'},
    {'name': 'var18', 'out_name': 'None', 'long_name': 'RESERVED FOR METVIEW. ', 'units': '-'},
    {'name': 'var19', 'out_name': 'None', 'long_name': 'RESERVED FOR METVIEW. ', 'units': '-'},
    {'name': 'var20', 'out_name': 'None', 'long_name': 'RESERVED FOR METVIEW. ', 'units': '-'},
    {'name': 'var21', 'out_name': 'UCTP', 'long_name': 'Unbalanced component of temperature ', 'units': 'K'},
    {'name': 'var22', 'out_name': 'UCLN', 'long_name': 'Unbalanced component of lnsp ', 'units': '-'},
    {'name': 'var23', 'out_name': 'UCDV', 'long_name': 'Unbalanced component of divergence', 'units': ' s**-1'},
    {'name': 'var24', 'out_name': 'None', 'long_name': 'Reserved for future unbalanced components', 'units': ' -'},
    {'name': 'var25', 'out_name': 'None', 'long_name': 'Reserved for future unbalanced components ', 'units': '-'},
    {'name': 'var26', 'out_name': 'CL', 'long_name': 'Lake cover ', 'units': '(0-1)'},
    {'name': 'var27', 'out_name': 'CVL', 'long_name': 'Low vegetation cover ', 'units': '(0-1)'},
    {'name': 'var28', 'out_name': 'CVH', 'long_name': 'High vegetation cover ', 'units': '(0-1)'},
    {'name': 'var29', 'out_name': 'TVL', 'long_name': 'Type of low vegetation ', 'units': '-'},
    {'name': 'var30', 'out_name': 'TVH', 'long_name': 'Type of high vegetation ', 'units': '-'},
    {'name': 'var31', 'out_name': 'CI', 'long_name': 'Sea ice cover ', 'units': '(0-1)'},
    {'name': 'var32', 'out_name': 'ASN', 'long_name': 'Snow albedo ', 'units': '(0-1)'},
    {'name': 'var33', 'out_name': 'RSN', 'long_name': 'Snow density ', 'units': 'kg**-3'},
    {'name': 'var34', 'out_name': 'SSTK', 'long_name': 'Sea surface temperature (absolute) ', 'units': 'K'},
    {'name': 'var35', 'out_name': 'ISTL1', 'long_name': 'Ice surface temperature layer 1 ', 'units': 'K'},
    {'name': 'var36', 'out_name': 'ISTL2', 'long_name': 'Ice surface temperature layer 2 ', 'units': 'K'},
    {'name': 'var37', 'out_name': 'ISTL3', 'long_name': 'Ice surface temperature layer 3 ', 'units': 'K'},
    {'name': 'var38', 'out_name': 'ISTL4', 'long_name': 'Ice surface temperature layer 4 ', 'units': 'K'},
    {'name': 'var39', 'out_name': 'SWVL1', 'long_name': 'Volumetric soil water layer 1 ', 'units': 'm**3 m**-3'},
    {'name': 'var40', 'out_name': 'SWVL2', 'long_name': 'Volumetric soil water layer 2 ', 'units': 'm**3 m**-3'},
    {'name': 'var41', 'out_name': 'SWVL3', 'long_name': 'Volumetric soil water layer 3 ', 'units': 'm**3 m**-3'},
    {'name': 'var42', 'out_name': 'SWVL4', 'long_name': 'Volumetric soil water layer 4 ', 'units': 'm**3 m**-3'},
    {'name': 'var43', 'out_name': 'SLT', 'long_name': 'Soil Type ', 'units': '-'},
    {'name': 'var44', 'out_name': '-', 'long_name': 'Unused ', 'units': '-'},
    {'name': 'var45', 'out_name': '-', 'long_name': 'Unused ', 'units': '-'},
    {'name': 'var46', 'out_name': '-', 'long_name': 'Unused ', 'units': '-'},
    {'name': 'var47', 'out_name': '-', 'long_name': 'Unused ', 'units': '-'},
    {'name': 'var48', 'out_name': '-', 'long_name': 'Unused ', 'units': '-'},
    {'name': 'var49', 'out_name': '-', 'long_name': 'Unused ', 'units': '-'},
    {'name': 'var50', 'out_name': '-', 'long_name': 'Unused ', 'units': '-'},
    {'name': 'var51', 'out_name': '-', 'long_name': 'Unused ', 'units': '-'},
    {'name': 'var52', 'out_name': '-', 'long_name': 'Unused ', 'units': '-'},
    {'name': 'var53', 'out_name': '-', 'long_name': 'Unused ', 'units': '-'},
    {'name': 'var54', 'out_name': '-', 'long_name': 'Unused ', 'units': '-'},
    {'name': 'var55', 'out_name': '-', 'long_name': 'Unused ', 'units': '-'},
    {'name': 'var56', 'out_name': '-', 'long_name': 'Unused ', 'units': '-'},
    {'name': 'var57', 'out_name': '-', 'long_name': 'Unused ', 'units': '-'},
    {'name': 'var58', 'out_name': '-', 'long_name': 'Unused ', 'units': '-'},
    {'name': 'var59', 'out_name': '-', 'long_name': 'Unused ', 'units': '-'},
    {'name': 'var60', 'out_name': '-', 'long_name': 'Unused ', 'units': '-'},
    {'name': 'var61', 'out_name': '-', 'long_name': 'Unused ', 'units': '-'},
    {'name': 'var62', 'out_name': '-', 'long_name': 'Unused ', 'units': '-'},
    {'name': 'var63', 'out_name': '-', 'long_name': 'Unused ', 'units': '-'},
    {'name': 'var64', 'out_name': '-', 'long_name': 'Unused ', 'units': '-'},
    {'name': 'var65', 'out_name': '-', 'long_name': 'Unused ', 'units': '-'},
    {'name': 'var66', 'out_name': '-', 'long_name': 'Unused ', 'units': '-'},
    {'name': 'var67', 'out_name': '-', 'long_name': 'Unused ', 'units': '-'},
    {'name': 'var68', 'out_name': '-', 'long_name': 'Unused ', 'units': '-'},
    {'name': 'var69', 'out_name': '-', 'long_name': 'Unused ', 'units': '-'},
    {'name': 'var70', 'out_name': '-', 'long_name': 'Unused ', 'units': '-'},
    {'name': 'var71', 'out_name': '-', 'long_name': 'Unused ', 'units': '-'},
    {'name': 'var72', 'out_name': '-', 'long_name': 'Unused ', 'units': '-'},
    {'name': 'var73', 'out_name': '-', 'long_name': 'Unused ', 'units': '-'},
    {'name': 'var74', 'out_name': '-', 'long_name': 'Unused ', 'units': '-'},
    {'name': 'var75', 'out_name': '-', 'long_name': 'Unused ', 'units': '-'},
    {'name': 'var76', 'out_name': '-', 'long_name': 'Unused ', 'units': '-'},
    {'name': 'var77', 'out_name': '-', 'long_name': 'Unused ', 'units': '-'},
    {'name': 'var78', 'out_name': '-', 'long_name': 'Unused ', 'units': '-'},
    {'name': 'var79', 'out_name': '-', 'long_name': 'Unused ', 'units': '-'},
    {'name': 'var80', 'out_name': '-', 'long_name': 'Unused ', 'units': '-'},
    {'name': 'var81', 'out_name': '-', 'long_name': 'Unused ', 'units': '-'},
    {'name': 'var82', 'out_name': '-', 'long_name': 'Unused ', 'units': '-'},
    {'name': 'var83', 'out_name': '-', 'long_name': 'Unused ', 'units': '-'},
    {'name': 'var84', 'out_name': '-', 'long_name': 'Unused ', 'units': '-'},
    {'name': 'var85', 'out_name': '-', 'long_name': 'Unused ', 'units': '-'},
    {'name': 'var86', 'out_name': '-', 'long_name': 'Unused ', 'units': '-'},
    {'name': 'var87', 'out_name': '-', 'long_name': 'Unused ', 'units': '-'},
    {'name': 'var88', 'out_name': '-', 'long_name': 'Unused ', 'units': '-'},
    {'name': 'var89', 'out_name': '-', 'long_name': 'Unused ', 'units': '-'},
    {'name': 'var90', 'out_name': '-', 'long_name': 'Unused ', 'units': '-'},
    {'name': 'var91', 'out_name': '-', 'long_name': 'Unused ', 'units': '-'},
    {'name': 'var92', 'out_name': '-', 'long_name': 'Unused ', 'units': '-'},
    {'name': 'var93', 'out_name': '-', 'long_name': 'Unused ', 'units': '-'},
    {'name': 'var94', 'out_name': '-', 'long_name': 'Unused ', 'units': '-'},
    {'name': 'var95', 'out_name': '-', 'long_name': 'Unused ', 'units': '-'},
    {'name': 'var96', 'out_name': '-', 'long_name': 'Unused ', 'units': '-'},
    {'name': 'var97', 'out_name': '-', 'long_name': 'Unused ', 'units': '-'},
    {'name': 'var98', 'out_name': '-', 'long_name': 'Unused ', 'units': '-'},
    {'name': 'var99', 'out_name': '-', 'long_name': 'Unused ', 'units': '-'},
    {'name': 'var100', 'out_name': '-', 'long_name': 'Unused ', 'units': '-'},
    {'name': 'var101', 'out_name': '-', 'long_name': 'Unused ', 'units': '-'},
    {'name': 'var102', 'out_name': '-', 'long_name': 'Unused ', 'units': '-'},
    {'name': 'var103', 'out_name': '-', 'long_name': 'Unused ', 'units': '-'},
    {'name': 'var104', 'out_name': '-', 'long_name': 'Unused ', 'units': '-'},
    {'name': 'var105', 'out_name': '-', 'long_name': 'Unused ', 'units': '-'},
    {'name': 'var106', 'out_name': '-', 'long_name': 'Unused ', 'units': '-'},
    {'name': 'var107', 'out_name': '-', 'long_name': 'Unused ', 'units': '-'},
    {'name': 'var108', 'out_name': '-', 'long_name': 'Unused ', 'units': '-'},
    {'name': 'var109', 'out_name': '-', 'long_name': 'Unused ', 'units': '-'},
    {'name': 'var110', 'out_name': '-', 'long_name': 'Unused ', 'units': '-'},
    {'name': 'var111', 'out_name': '-', 'long_name': 'Unused ', 'units': '-'},
    {'name': 'var112', 'out_name': '-', 'long_name': 'Unused ', 'units': '-'},
    {'name': 'var113', 'out_name': '-', 'long_name': 'Unused ', 'units': '-'},
    {'name': 'var114', 'out_name': '-', 'long_name': 'Unused ', 'units': '-'},
    {'name': 'var115', 'out_name': '-', 'long_name': 'Unused ', 'units': '-'},
    {'name': 'var116', 'out_name': '-', 'long_name': 'Unused ', 'units': '-'},
    {'name': 'var117', 'out_name': '-', 'long_name': 'Unused ', 'units': '-'},
    {'name': 'var118', 'out_name': '-', 'long_name': 'Unused ', 'units': '-'},
    {'name': 'var119', 'out_name': '-', 'long_name': 'Unused ', 'units': '-'},
    {'name': 'var120', 'out_name': '-', 'long_name': 'Unused ', 'units': '-'},
    {'name': 'var121', 'out_name': '-', 'long_name': 'Unused ', 'units': '-'},
    {'name': 'var122', 'out_name': '-', 'long_name': 'Unused ', 'units': '-'},
    {'name': 'var123', 'out_name': '-', 'long_name': 'Unused ', 'units': '-'},
    {'name': 'var124', 'out_name': '-', 'long_name': 'Unused ', 'units': '-'},
    {'name': 'var125', 'out_name': '-', 'long_name': 'Unused ', 'units': '-'},
    {'name': 'var126', 'out_name': '-', 'long_name': 'Unused ', 'units': '-'},
    {'name': 'var127', 'out_name': 'AT', 'long_name': 'Atmospheric tide ', 'units': '-'},
    {'name': 'var128', 'out_name': 'BV', 'long_name': 'Budget values ', 'units': '-'},
    {'name': 'var129', 'out_name': 'Z', 'long_name': 'Geopotential (at the surface = orography) ',
     'units': 'm**2 s**-2'},
    {'name': 'var130', 'out_name': 'T', 'long_name': 'Temperature ', 'units': 'K'},
    {'name': 'var131', 'out_name': 'U', 'long_name': 'U-velocity ', 'units': 'm s**-1'},
    {'name': 'var132', 'out_name': 'V', 'long_name': 'V-velocity ', 'units': 'm s**-1'},
    {'name': 'var133', 'out_name': 'Q', 'long_name': 'Specific humidity ', 'units': 'kg kg**-1'},
    {'name': 'var134', 'out_name': 'SP', 'long_name': 'Surface pressure ', 'units': 'Pa'},
    {'name': 'var135', 'out_name': 'W', 'long_name': 'Vertical velocity ', 'units': 'Pa s**-1'},
    {'name': 'var136', 'out_name': '-', 'long_name': 'Unused ', 'units': '-'},
    {'name': 'var137', 'out_name': 'PWC', 'long_name': 'Precipitable water content ', 'units': 'kg m**-2'},
    {'name': 'var138', 'out_name': 'VO', 'long_name': 'Vorticity (relative) ', 'units': 's**-1'},
    {'name': 'var139', 'out_name': 'STL1', 'long_name': 'Soil temperature level 1 ', 'units': 'K'},
    {'name': 'var140', 'out_name': 'SWL1', 'long_name': 'Soil wetness level 1 ', 'units': 'm'},
    {'name': 'var141', 'out_name': 'SD', 'long_name': 'Snow depth ', 'units': 'm (of water)'},
    {'name': 'var142', 'out_name': 'LSP', 'long_name': 'Large scale precipitation ', 'units': 'kg m**-2 s**-1'},
    {'name': 'var143', 'out_name': 'CP', 'long_name': 'Convective precipitation ', 'units': 'kg m**-2 s**-1'},
    {'name': 'var144', 'out_name': 'SF', 'long_name': 'Snow fall ', 'units': 'kg m**-2 s**-1'},
    {'name': 'var145', 'out_name': 'BLD', 'long_name': 'Boundary layer dissipation ', 'units': 'W m**-2 s'},
    {'name': 'var146', 'out_name': 'SSHF', 'long_name': 'Surface sensible heat flux ', 'units': 'W m**-2 s'},
    {'name': 'var147', 'out_name': 'SLHF', 'long_name': 'Surface latent heat flux ', 'units': 'W m**-2 s'},
    {'name': 'var148', 'out_name': '-', 'long_name': 'Unused ', 'units': '-'},
    {'name': 'var149', 'out_name': '-', 'long_name': 'Unused ', 'units': '-'},
    {'name': 'var150', 'out_name': '-', 'long_name': 'Unused ', 'units': '-'},
    {'name': 'var151', 'out_name': 'MSL', 'long_name': 'Mean sea level pressure', 'units': 'Pa'},
    {'name': 'var152', 'out_name': 'LNSP', 'long_name': 'Ln surface pressure', 'units': ' -'},
    {'name': 'var153', 'out_name': '-', 'long_name': 'Unused ', 'units': '-'},
    {'name': 'var154', 'out_name': '-', 'long_name': 'Unused ', 'units': '-'},
    {'name': 'var155', 'out_name': 'D', 'long_name': 'Divergence', 'units': 's**-1'},
    {'name': 'var156', 'out_name': 'GH', 'long_name': 'Height (geopotential)', 'units': 'm'},
    {'name': 'var157', 'out_name': 'R', 'long_name': 'Relative humidity', 'units': '(0 - 1)'},
    {'name': 'var158', 'out_name': 'TSP', 'long_name': 'Tendency of surface pressure ', 'units': 'Pa s**-1'},
    {'name': 'var159', 'out_name': '-', 'long_name': 'Unused ', 'units': '-'},
    {'name': 'var160', 'out_name': '-', 'long_name': 'Unused ', 'units': '-'},
    {'name': 'var161', 'out_name': '-', 'long_name': 'Unused ', 'units': '-'},
    {'name': 'var162', 'out_name': '-', 'long_name': 'Unused ', 'units': '-'},
    {'name': 'var163', 'out_name': '-', 'long_name': 'Unused ', 'units': '-'},
    {'name': 'var164', 'out_name': 'TCC', 'long_name': 'Total cloud cover ', 'units': '(0 - 1)'},
    {'name': 'var165', 'out_name': '10U', 'long_name': '10 metre u wind component', 'units': ' m s**-1'},
    {'name': 'var166', 'out_name': '10V', 'long_name': '10 metre v wind component ', 'units': 'm s**-1'},
    {'name': 'var167', 'out_name': '2T', 'long_name': '2 metre temperature ', 'units': 'K'},
    {'name': 'var168', 'out_name': '2D', 'long_name': '2 metre dewpoint temperature', 'units': ' K'},
    {'name': 'var169', 'out_name': '-', 'long_name': 'Unused ', 'units': '-'},
    {'name': 'var170', 'out_name': 'STL2', 'long_name': 'Soil temperature level 2 ', 'units': 'K'},
    {'name': 'var171', 'out_name': 'SWL2', 'long_name': 'Soil wetness level 2 ', 'units': 'm'},
    {'name': 'var172', 'out_name': 'LSM', 'long_name': 'Land/sea mask ', 'units': '(0 - 1)'},
    {'name': 'var173', 'out_name': 'SR', 'long_name': 'Surface roughness ', 'units': 'm'},
    {'name': 'var174', 'out_name': 'AL', 'long_name': 'Albedo ', 'units': '(0 - 1)'},
    {'name': 'var175', 'out_name': '-', 'long_name': 'Unused ', 'units': '-'},
    {'name': 'var176', 'out_name': 'SSR', 'long_name': 'Surface solar radiation ', 'units': 'W m**-2'},
    {'name': 'var177', 'out_name': 'STR', 'long_name': 'Surface thermal radiation ', 'units': 'W m**-2'},
    {'name': 'var178', 'out_name': 'TSR', 'long_name': 'Top solar radiation ', 'units': 'W m**-2'},
    {'name': 'var179', 'out_name': 'TTR', 'long_name': 'Top thermal radiation ', 'units': 'W m**-2'},
    {'name': 'var180', 'out_name': 'EWSS', 'long_name': 'East/west surface stress', 'units': ' N m**-2 s**-1'},
    {'name': 'var181', 'out_name': 'NSSS', 'long_name': 'North/south surface stress ', 'units': 'N m**-2 s**-1'},
    {'name': 'var182', 'out_name': 'E', 'long_name': 'Evaporation ', 'units': 'kg m**-2 s**-1'},
    {'name': 'var183', 'out_name': 'STL3', 'long_name': 'Soil temperature level 3 ', 'units': 'K'},
    {'name': 'var184', 'out_name': 'SWL3', 'long_name': 'Soil wetness level 3 ', 'units': 'm'},
    {'name': 'var185', 'out_name': 'CCC', 'long_name': 'Convective cloud cover ', 'units': '(0 - 1)'},
    {'name': 'var186', 'out_name': 'LCC', 'long_name': 'Low cloud cover ', 'units': '(0 - 1)'},
    {'name': 'var187', 'out_name': 'MCC', 'long_name': 'Medium cloud cover ', 'units': '(0 - 1)'},
    {'name': 'var188', 'out_name': 'HCC', 'long_name': 'High cloud cover ', 'units': '(0 - 1)'},
    {'name': 'var189', 'out_name': '-', 'long_name': 'Unused ', 'units': '-'},
    {'name': 'var190', 'out_name': 'EWOV', 'long_name': 'EW component of sub-grid_routines scale orographic variance',
     'units': ' m**2'},
    {'name': 'var191', 'out_name': 'NSOV', 'long_name': 'NS component of sub-grid_routines scale orographic variance',
     'units': ' m**2'},
    {'name': 'var192', 'out_name': 'NWOV', 'long_name': 'NWSE component sub-grid_routines scale orographic variance ',
     'units': 'm**2'},
    {'name': 'var193', 'out_name': 'NEOV', 'long_name': 'NESW component sub-grid_routines scale orographic variance',
     'units': ' m**2'},
    {'name': 'var194', 'out_name': '-', 'long_name': 'Unused ', 'units': '-'},
    {'name': 'var195', 'out_name': 'LGWS', 'long_name': 'Latitudinal component of gravity wave stress ',
     'units': 'N m**-2 s'},
    {'name': 'var196', 'out_name': 'MGWS', 'long_name': 'Meridional component of gravity wave stress ',
     'units': 'N m**-2 s'},
    {'name': 'var197', 'out_name': 'GWD', 'long_name': 'Gravity wave dissipation', 'units': ' W m**-2 s'},
    {'name': 'var198', 'out_name': 'SRC', 'long_name': 'Skin reservoir content ', 'units': 'm (of water)'},
    {'name': 'var199', 'out_name': 'VEG', 'long_name': 'Percentage of vegetation ', 'units': '%'},
    {'name': 'var200', 'out_name': 'VSO', 'long_name': 'Variance of sub-grid_routines scale orography ',
     'units': 'm**2'},
    {'name': 'var201', 'out_name': 'MX2T', 'long_name': 'Max temp.2m during averaging time ', 'units': 'K'},
    {'name': 'var202', 'out_name': 'MN2T', 'long_name': 'Min temp.2m during averaging time', 'units': ' K'},
    {'name': 'var203', 'out_name': '-', 'long_name': 'Unused ', 'units': '-'},
    {'name': 'var204', 'out_name': 'PAW', 'long_name': 'Precip. analysis weights ', 'units': '-'},
    {'name': 'var205', 'out_name': 'RO', 'long_name': 'Runoff ', 'units': 'kg m**-2 s**-1'},
    {'name': 'var206', 'out_name': 'ZZ', 'long_name': 'St.Dev. of Geopotential ', 'units': 'm**2 s**-2'},
    {'name': 'var207', 'out_name': 'TZ', 'long_name': 'Covar Temp & Geopotential', 'units': ' K m**2 s**-2'},
    {'name': 'var208', 'out_name': 'TT', 'long_name': 'St.Dev. of Temperature', 'units': ' K'},
    {'name': 'var209', 'out_name': 'QZ', 'long_name': 'Covar Sp.Hum. & Geopotential', 'units': ' m**2 s**-2'},
    {'name': 'var210', 'out_name': 'QT', 'long_name': 'Covar Sp.Hum & Temp. ', 'units': 'K'},
    {'name': 'var211', 'out_name': 'QQ', 'long_name': 'St.Dev. of Specific humidity ', 'units': '(0 - 1)'},
    {'name': 'var212', 'out_name': 'UZ', 'long_name': 'Covar U-comp. & Geopotential', 'units': ' m**3 s**-3'},
    {'name': 'var213', 'out_name': 'UT', 'long_name': 'Covar U-comp. & Temp.', 'units': ' K m s**-1'},
    {'name': 'var214', 'out_name': 'UQ', 'long_name': 'Covar U-comp. & Sp.Hum.', 'units': ' m s**-1'},
    {'name': 'var215', 'out_name': 'UU', 'long_name': 'St.Dev. of U-velocity ', 'units': 'm s**-1'},
    {'name': 'var216', 'out_name': 'VZ', 'long_name': 'Covar V-comp. & ', 'units': 'Geopotential m**3 s**-3'},
    {'name': 'var217', 'out_name': 'VT', 'long_name': 'Covar V-comp. & Temp.', 'units': ' K m s**-1'},
    {'name': 'var218', 'out_name': 'VQ', 'long_name': 'Covar V-comp. & Sp.Hum.', 'units': ' m s**-1'},
    {'name': 'var219', 'out_name': 'VU', 'long_name': 'Covar V-comp. & U-comp ', 'units': 'm**2 s**-2'},
    {'name': 'var220', 'out_name': 'VV', 'long_name': 'St.Dev. of V-comp', 'units': ' m s**-1'},
    {'name': 'var221', 'out_name': 'WZ', 'long_name': 'Covar W-comp. & Geopotential ', 'units': 'Pa m**2 s**-3'},
    {'name': 'var222', 'out_name': 'WT', 'long_name': 'Covar W-comp. & Temp.', 'units': ' K Pa s**-1'},
    {'name': 'var223', 'out_name': 'WQ', 'long_name': 'Covar W-comp. & Sp.Hum.', 'units': ' Pa s**-1'},
    {'name': 'var224', 'out_name': 'WU', 'long_name': 'Covar W-comp. & U-comp.', 'units': ' Pa m s**-2'},
    {'name': 'var225', 'out_name': 'WV', 'long_name': 'Covar W-comp. & V-comp.', 'units': ' Pa m s**-2'},
    {'name': 'var226', 'out_name': 'Something', 'long_name': 'Unused ', 'units': '-'},
    {'name': 'var227', 'out_name': 'WW', 'long_name': 'St.Dev. of Vertical velocity Pa', 'units': 's**-1'},
    {'name': 'var228', 'out_name': 'Something_else', 'long_name': 'Unused ', 'units': '-'},
    {'name': 'var229', 'out_name': 'TP', 'long_name': 'Total precipitation', 'units': 'm'},
    {'name': 'var230', 'out_name': 'IEWS', 'long_name': 'Instantaneous X surface stress', 'units': 'N m**-2'},
    {'name': 'var231', 'out_name': 'INSS', 'long_name': 'Instantaneous Y surface stress', 'units': 'N m**-2'},
    {'name': 'var232', 'out_name': 'ISHF', 'long_name': 'Instantaneous surface Heat Flux', 'units': 'W m**-2'},
    {'name': 'var233', 'out_name': 'IE', 'long_name': 'Instantaneous Moisture Flux (evaporation)',
     'units': 'kg m**-2 s**-1'},
    {'name': 'var234', 'out_name': 'ASQ', 'long_name': 'Apparent Surface Humidity ', 'units': 'kg kg**-1'},
    {'name': 'var235', 'out_name': 'LSRH', 'long_name': 'Logarithm of surface roughness length for heat',
     'units': '-'},
    {'name': 'var236', 'out_name': 'SKT', 'long_name': 'Skin Temperature', 'units': 'K'},
    {'name': 'var237', 'out_name': 'STL4', 'long_name': 'Soil temperature level 4 ', 'units': 'K'},
    {'name': 'var238', 'out_name': 'SWL4', 'long_name': 'Soil wetness level 4 ', 'units': 'm'},
    {'name': 'var239', 'out_name': 'TSN', 'long_name': 'Temperature of snow layer ', 'units': 'K'},
    {'name': 'var240', 'out_name': 'CSF', 'long_name': 'Convective snow-fall ', 'units': 'kg m**-2 s**-1'},
    {'name': 'var241', 'out_name': 'LSF', 'long_name': 'Large scale snow-fall ', 'units': 'kg m**-2 s**-1'},
    {'name': 'var242', 'out_name': 'CLWC', 'long_name': 'Cloud liquid water content ', 'units': 'kg kg**-1'},
    {'name': 'var243', 'out_name': 'CC', 'long_name': 'Cloud cover (at given level) ', 'units': '(0 - 1)'},
    {'name': 'var244', 'out_name': 'FAL', 'long_name': 'Forecast albedo ', 'units': '-'},
    {'name': 'var245', 'out_name': 'FSR', 'long_name': 'Forecast surface roughness ', 'units': 'm'},
    {'name': 'var246', 'out_name': 'FLSR', 'long_name': 'Forecast logarithm of surface roughness for heat. ',
     'units': '-'},
    {'name': 'var247', 'out_name': '10WS', 'long_name': '10m. Windspeed (irresp of dir.) ', 'units': 'm s**-1'},
    {'name': 'var248', 'out_name': 'MOFL', 'long_name': 'Momentum flux (irresp of dir.) ', 'units': 'N m**-2'},
    {'name': 'var249', 'out_name': 'HSD', 'long_name': 'Heaviside ("beta") function', 'units': '(0 - 1)'},
    {'name': 'var250', 'out_name': '-', 'long_name': 'Unused ', 'units': '-'},
    {'name': 'var251', 'out_name': '-', 'long_name': 'Unused ', 'units': '-'},
    {'name': 'var252', 'out_name': '-', 'long_name': 'Unused ', 'units': '-'},
    {'name': 'var253', 'out_name': '-', 'long_name': 'Unused ', 'units': '-'},
    {'name': 'var254', 'out_name': '-', 'long_name': 'Unused ', 'units': '-'},
    {'name': 'var255', 'out_name': '-', 'long_name': 'Unused ', 'units': '-'}]


def get_para_table_file_var(para_table: dict, var_name: str, source: str, index: int) -> str:
    """
    Build a CDO parameter table of a single variable from a dict
    :param para_table: dict which contains conversion information
    must include 'name' the name in the input file
    must include 'out_name' the name in the output file
    :param var_name: variable which is desired. is an 'out_name' in para_table
    :param source: used to add attribute
    :param index: used to add attribute of source
    :return: string of the parameter table
    """
    assert source[-4:] == '.grb'
    try:
        var_dict = next(x for x in para_table if x['out_name'] == var_name)
        parameter_desc = "\n&parameter\n" + \
                         '\n'.join([key + ' = ' + value for key, value in var_dict.items()]) + \
                         '\n' + var_dict['out_name'] + 'file' + str(index) + " = " + source + '\n/\n'
        filename = 'para_table_' + '_'.join(source[:-4].split('/')) + '_' + var_name + '.txt'
        text_file = open(filename, 'w')
        text_file.write(parameter_desc)
        text_file.close()
        return filename
    except StopIteration:
        sys.exit('ERROR: Varname not in parameter table!')


def get_para_table_file_vars(para_table: dict, var_list: list, file: str, all_files: str) -> str:
    """
    Build a CDO parameter table of a single variable from a dict
    :param para_table: dict which contains conversion information
    must include 'name' the name in the input file
    must include 'out_name' the name in the output file
    :param var_list: list of variables to extract from the file, paired if the name in the file is different
    :param file: used to add attribute
    :param all_files: list of all files used as a source for these variables
    :return: string of the parameter table
    """
    assert file[-4:] == '.grb'
    filename = 'para_table_' + '_'.join(file[:-4].split('/')) + '.txt'
    text_file = open(filename, 'w')
    var_dicts = [x for x in para_table if x['out_name'] in var_list]
    for var_dict in var_dicts:
        try:
            parameter_desc = "\n&parameter\n" + \
                             '\n'.join([key + ' = ' + value for key, value in var_dict.items()]) + \
                             '\n' + var_dict['out_name'] + 'files = ' + all_files + '\n/\n'
            text_file.write(parameter_desc)
        except StopIteration:
            text_file.close()
            sys.exit('ERROR: Varname not in parameter table!')
    text_file.close()
    return filename


def get_para_table_file(para_table: dict, source: str) -> str:
    """
    Build a CDO parameter table to be used with setpartabn from a dictionary
    :param para_table: dict which contains conversion information
    must include 'name' the name in the input file
    must include 'out_name' the name in the output file
    :param source: Used to add an attribute detailing the source file
    :return: string of the parameter table
    """
    assert source[-4:] == '.grb'
    found_names = []
    parameter_desc = ""
    for a_dict in para_table:
        if a_dict['out_name'] == '-' or a_dict['out_name'] in found_names:
            continue
        parameter_desc += "\n&parameter\n" + \
                          '\n'.join([key + ' = ' + value for key, value in a_dict.items()]) + \
                          '\n' + a_dict['out_name'] + 'source = ' + source + '\n/\n'
        found_names += [a_dict['out_name']]
    filename = 'para_table_' + '_'.join(source[:-4].split('/')) + '.txt'
    text_file = open(filename, 'w')
    text_file.write(parameter_desc)
    text_file.close()
    return filename
