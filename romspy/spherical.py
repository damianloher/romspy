import numpy as np

def spheredist(lon,lat):
    """
    Gives the distance in kilometers between
    successive points in the vectors lon and lat, computed
    using the Haversine formula on a spherical earth of radius 6371.009km. 
    :param lon: numpy 1D array of longitudes (deg)
    :param lat: numpy 1D array of latitudes (deg)
    return: numpy 1D array of distances between corresponding points
       in km
    """
    if isinstance(lon, list):
        lon = np.array(lon)
    if isinstance(lat, list):
        lat = np.array(lat)
    if lon.ndim != 1 or lat.ndim != 1:
        msg = 'lon and lat need to be 1D numpy arrays'
        raise ValueError(msg)
    if lon.shape[0] != lat.shape[0]:
        msg = 'lon and lat need to be of the same length'
        raise ValueError(msg)
    if lon.shape[0] < 2:
        msg = 'at least 2 points need to be given'
        raise ValueError(msg)
    deg2rad = np.pi/180
    # Mean equatorial radius according to IUGG [km]:
    er = 6371.009
    lon = lon*deg2rad
    lat = lat*deg2rad
    # Starting latitudes:
    lat1 = lat[:-1]
    # Ending latitudes:
    lat2 = lat[1:]
    dlon = np.diff(lon)
    dlat = np.diff(lat)
    # Compute distances in km:
    dist = 2*er*np.arcsin(np.sqrt( np.sin(0.5*dlat)**2 + 
            np.cos(lat1)*np.cos(lat2)*np.sin(0.5*dlon)**2 ))
    return dist

