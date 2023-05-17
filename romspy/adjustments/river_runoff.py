import netCDF4
import os
import numpy as np
from romspy import UP_data_paths
from romspy import spherical
import sys

"""
Freshwater flux correction due to river runoff. Based on make_river_freshwater.m
from the Romstools.
"""

def make_river_freshwater(grd_file,roms_setup,out_dir,**kwargs):
    """
    River runoff on the ROMS grid, spread over a specified number of
    ROMS grid cells.
    :param grd_file: path to ROMS grid file
    :param out_dir: output folder
    Can have any of the following optional arguments:
        verbose - whether to print runtime information - default false
    """
    # Check if output file exists: assume nothing needs to be done
    # if this is the case:
    r_tag = roms_setup.lower()
    out_file = '{}/{}_1day_river_frc.nc'.format(out_dir,r_tag)
    if os.path.exists(out_file):
        print('found river runoff on ROMS grid: '+out_file)
        return out_file
    verbose = kwargs.get('verbose', False)
    spreading = kwargs.get('spreading', 'gaussian')
    if verbose:
        print('compute river runoff on ROMS grid')
    from shapely.geometry import Point
    from shapely.geometry.polygon import Polygon
    import scipy
    # Load Dai et al. (2009) data:
    infile = UP_data_paths.Dai_river_runoff_dir + 'Dai2009_river_discharge_monthly_clim.nc'
    nc = netCDF4.Dataset(infile,'r')
    Lon_dat = nc.variables['lon_mou'][:]
    Lat_dat = nc.variables['lat_mou'][:]
    flow = nc.variables['FLOW_CLIM'][:]  # shape: time_clim, station
    nc.close()
    # Load variables related to the ROMS grid:
    nc = netCDF4.Dataset(grd_file,'r')
    lon_rho = nc.variables['lon_rho'][:]
    lat_rho = nc.variables['lat_rho'][:]
    mask_rho = nc.variables['mask_rho'][:]
    pm = nc.variables['pm'][:]
    pn = nc.variables['pn'][:]
    nc.close()
    area = 1/(pm*pn)  # m^2
    ny, nx = mask_rho.shape
    nt = flow.shape[0]
    nriv = flow.shape[1]
    river_in = []
    # Figure out which river mouthes are in the ROMS domain and
    # do some manual adjustments based on the ROMS setup:
    if roms_setup.lower() == 'pactcs30':
        # Make sure Lon_dat is in [0,360):
        Lon_dat[Lon_dat<0] += 360
        # Construct polygon representing the boundary of the
        # ROMS domain:
        plist = []
        for i in range(nx):
            plist.append((lon_rho[0,i],lat_rho[0,i]))
        for j in range(1,ny):
            plist.append((lon_rho[j,-1],lat_rho[j,-1]))
        for i in range(nx-2,-1,-1):
            plist.append((lon_rho[-1,i],lat_rho[-1,i]))
        for j in range(ny-2,0,-1):
            plist.append((lon_rho[j,0],lat_rho[j,0]))
        polygon = Polygon(plist)
        if verbose:
            print('   Figure out which rivers are in the ROMS domain...')
            sys.stdout.flush()
        for r in range(nriv):
            if polygon.contains(Point(Lon_dat[r],Lat_dat[r])):
                riv_in = True
                # Exclude some rivers:
                if Lat_dat[r]>65:
                    # North Pacific:
                    riv_in = False
                elif Lon_dat[r]>300:
                    riv_in = False
                elif (Lon_dat[r]>260 and Lon_dat[r]<270) and (Lat_dat[r]>17 and Lat_dat[r]<23):
                    # Central America
                    riv_in = False
                elif (Lon_dat[r]>280 and Lon_dat[r]<282) and (Lat_dat[r]>9 and Lat_dat[r]<11):
                    # Central America
                    riv_on = False
            else:
                riv_in = False
                if (Lon_dat[r]>236 and Lon_dat[r]<240) and (Lat_dat[r]>37 and Lat_dat[r]<42):
                    # NW Pacific:
                    riv_in = True
            river_in.append(riv_in)
    else:
        msg = 'river correction is not implemented for this ROMS setup: {}'.format(roms_setup)
        raise NotImplementedError(msg)
    # Go through river data and find the cell on the ROMS grid that is
    # closest to the river mouth:
    score = np.zeros((ny,nx), dtype=np.float32)
    flow_ROMSgrid = np.zeros((nt,ny,nx), dtype=np.float32)
    if verbose:
        print('   Find nearest ROMS grid cell for each river...')
        sys.stdout.flush()
    for r in range(nriv):
        if verbose:
            print('     r = {}/{}'.format(r,nriv))
            sys.stdout.flush()
        if river_in[r]:
            for j in range(ny):
                for i in range(nx):
                    # Calculate a score/distance for each ocean point:
                    if mask_rho[j,i] == 0:
                        score[j,i] = 1e+23
                    else:
                        score[j,i] = np.sqrt((Lon_dat[r]-lon_rho[j,i])**2 + (Lat_dat[r]-lat_rho[j,i])**2)
                    # Find cell with minimum distance:
                    iy, ix = np.unravel_index(np.argmin(score, axis=None), score.shape)
                    for t in range(nt):
                        if not np.isnan(flow[t,r]):
                            if r == 17 and r_tag.lower() == 'pactcs30':
                                # Make sure Amur is added in Sea of Okhotsk:
                                flow_ROMSgrid[t,105,88] += flow[t,r]
                            else:
                                flow_ROMSgrid[t,iy,ix] += flow[t,r]
    # Spread river runoff over a certain area of grid points depending on the
    # magnitude of the freshwater discharge in either linear or gaussian
    # fashion:
    min_extraboxes = 9
    if verbose:
        print('   Spreading river discharge on area, weighted by distance squared...')
        sys.stdout.flush()
    FLOW_spread = np.zeros((nt,ny,nx), dtype=np.float32)
    count = 0
    countextrab = 0
    check = False
    for t in range(nt):
        for j in range(ny):
            for i in range(nx):
                if flow_ROMSgrid[t,j,i] > 0:
                    FLOW_check = 0
                    # Number of surrounding indices to take: min=3, max=7 depending on absolute water runoff  
                    # (the max option was for the case of the "constant radius" spread)
                    tmp = int(np.floor(np.log10(flow_ROMSgrid[t,j,i])))
                    extraboxes = tmp + min_extraboxes
                    if verbose and (tmp < 0):
                        # Make sure that you have a decent min number of extraboxes
                        print('     ** warning: log10(FLOW)<0 ** --> log10(FLOW) = {} --> extraboxes = {}'.\
                            format(tmp,extraboxes))
                    # Define box indices:
                    istart = max(0, i-extraboxes)
                    jstart = max(0, j-extraboxes)
                    iend = min(nx-1, i+extraboxes)
                    jend = min(ny-1, j+extraboxes)
                    dist = np.nan*np.ones((jend-jstart+1,iend-istart+1), dtype=np.float32)
                    # Computes distances from cell (i,j):
                    for jbox in range(jstart,jend+1):
                        for ibox in range(istart,iend+1):
                            if mask_rho[jbox,ibox] == 1:
                                dist[jbox-jstart,ibox-istart] = spherical.spheredist([lon_rho[jbox,ibox], lon_rho[j,i]], \
                                    [lat_rho[jbox,ibox], lat_rho[j,i]])[0]
                    # Construct weighting array:
                    maxdist = dist.max()
                    if spreading == 'gaussian':
                        # Weighting with the square of the distance:
                        distweighted = np.exp(-(dist/(maxdist/2))**2)
                    else:
                        distweighted = 1-(dist/maxdist)
                    weights = distweighted/np.nansum(distweighted)
                    # Spreads the weighted data across the box:
                    for jbox in range(jstart,jend+1):
                        for ibox in range(istart,iend+1):
                            if mask_rho[jbox,ibox] == 1:
                                FLOW_spread[t,jbox,ibox] = FLOW_spread[t,jbox,ibox] + \
                                    flow_ROMSgrid[t,j,i]*weights[jbox-jstart,ibox-istart]
                                if check:
                                    FLOW_check = FLOW_check + \
                                        flow_ROMSgrid[t,j,i]*weights[jbox-jstart,ibox-istart]
                    if check:
                        print('FLOW_check/flow_ROMSgrid[t,j,i] = {}'.format(FLOW_check/flow_ROMSgrid[t,j,i]))
    # Convert river freshwater units (m3/s to cm/day)
    # (correct direction: neg. swflux for added freshwater!!!)
    for t in range(nt):
        FLOW_spread[t,:] = -(FLOW_spread[t,:]*86400/area) * 100
    #FLOW_spread = -((FLOW_spread.*86400)./repmat(gridarea,[1,1,t])).*100
    # Interpolate data from monthly to daily:
    molen = np.array([31,29,31,30,31,30,31,31,30,31,30,31])
    molen2 = np.zeros((13))
    molen2[1:] = molen
    midmo = np.cumsum(molen2[:-1])+0.5*molen
    # Extend midmo to prevent extrapolation:
    midmo = np.concatenate(([-15.5], midmo, [381.5]))
    mid_daily = np.arange(0.5,365.5)
    tmp = np.zeros((len(mid_daily),ny,nx), dtype=np.float32)
    if verbose:
        print('   compute daily river runoff from monthly...')
        sys.stdout.flush()
    for j in range(ny):
        for i in range(nx):
            tmp_val = np.concatenate(([FLOW_spread[-1,j,i]], FLOW_spread[:,j,i], [FLOW_spread[0,j,i]]))
            func = scipy.interpolate.interp1d(midmo, tmp_val)
            tmp[:,j,i] = func(mid_daily)

    # Write data to output file:
    if verbose:
        print('   write data to river input file: {}'.format(out_file))
    nc_out = netCDF4.Dataset(out_file,'w')
    nc_out.createDimension('xi_rho', size=nx)
    nc_out.createDimension('eta_rho', size=ny)
    nc_out.createDimension('swf_time', size=len(mid_daily))
    # Create variables:
    vobj = nc_out.createVariable('swf_time', 'f4', ('swf_time',))
    vobj.cycle_length = 365.0
    vobj.calendar = "365_day"
    vobj[:] = mid_daily
    vobj = nc_out.createVariable('swflux', 'f4', ('swf_time','eta_rho','xi_rho',))
    vobj.long_name = "Freshwater flux from river runoff"
    vobj.description = "Dai et al. (2009) climatological runoff data"
    vobj.units = "cm day-1"
    vobj[:] = tmp
    nc_out.close()
    # Return path to output file:
    return out_file