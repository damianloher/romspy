#include <print>
#include <vector>
#include <cmath>
#include <iostream>
#include <netcdf>
#include <vector>
#include "nanoflann.hpp"
#include "roms_forcing_utils.h"

using namespace nanoflann;
using namespace netCDF;

constexpr double DEG2RAD = M_PI / 180.0;

// =====================================================
// Point cloud in 3D Cartesian coordinates
// =====================================================
struct PointCloud
{
    std::vector<double> x, y, z, val;

    inline size_t kdtree_get_point_count() const
    {
        return x.size();
    }

    inline double kdtree_get_pt(const size_t idx, int dim) const
    {
        if (dim == 0) return x[idx];
        if (dim == 1) return y[idx];
        return z[idx];
    }

    template <class BBOX>
    bool kdtree_get_bbox(BBOX&) const
    {
        return false;
    }
};

// =====================================================
// Earth interpolator
// =====================================================
class EarthInterpolant
{
public:
    using KDTree = KDTreeSingleIndexAdaptor<
        L2_Simple_Adaptor<double, PointCloud>,
        PointCloud,
        3
    >;

    PointCloud pts;
    KDTree index;

    size_t k;
    double power;

    EarthInterpolant(
        const std::vector<double>& lat_deg,
        const std::vector<double>& lon_deg,
        const std::vector<double>& values,
        size_t k_ = 6,
        double power_ = 2.0)
        :
        index(3, pts, KDTreeSingleIndexAdaptorParams(10)),
        k(k_),
        power(power_)
    {
        const size_t N = lat_deg.size();

        pts.x.resize(N);
        pts.y.resize(N);
        pts.z.resize(N);
        pts.val = values;

        for (size_t i = 0; i < N; ++i)
        {
            double lat = lat_deg[i] * DEG2RAD;
            double lon = lon_deg[i] * DEG2RAD;

            pts.x[i] = std::cos(lat) * std::cos(lon);
            pts.y[i] = std::cos(lat) * std::sin(lon);
            pts.z[i] = std::sin(lat);
        }

        index.buildIndex();
    }

    double operator()(double latq_deg, double lonq_deg) const
    {
        double lat = latq_deg * DEG2RAD;
        double lon = lonq_deg * DEG2RAD;

        double query_pt[3] = {
            std::cos(lat) * std::cos(lon),
            std::cos(lat) * std::sin(lon),
            std::sin(lat)
        };

        std::vector<uint32_t> ret_index(k);
        std::vector<double> out_dist_sqr(k);

        KNNResultSet<double, uint32_t> resultSet(k);
        resultSet.init(ret_index.data(), out_dist_sqr.data());

        index.findNeighbors(resultSet, query_pt);

        double wsum = 0.0;
        double zsum = 0.0;

        for (size_t i = 0; i < k; ++i)
        {
            double d = std::sqrt(out_dist_sqr[i]);
            d = std::max(d, 1e-12);   // avoid division by zero

            double w = 1.0 / std::pow(d, power);

            wsum += w;
            zsum += w * pts.val[ret_index[i]];
        }

        return zsum / wsum;
    }
};

void ROMS_forcing::era_extrapolate_land(std::string ERA_file, std::vector<std::string> variables)
{
    std::println("\n------- ROMS_forcing::era_extrapolate_land -------");
    std::println("ERA file:  {}", ERA_file);
    NcFile nc;
    NcVar v, v_lon, v_lat;
    size_t nt, ny, nx, nxy;
    // Open Netcdf file:
    nc.open(ERA_file, NcFile::FileMode::write);
    nt = nc.getDim("time").getSize();
    v_lon = nc.getVar("lon");
    v_lat = nc.getVar("lat");
    for (const auto &var : variables) {
        std::cout << "var = " << var << '\n';
        std::cout.flush();
        v = nc.getVar(var);
        // Allocate buffers for data, latitudes and longitudes:
        nxy = nx*ny;
        std::vector<double> xy_slice(nxy);
        std::vector<double> lon(nx), lat(ny);
        std::vector<double> xy_slice_valid;
        std::vector<double> lon_valid, lat_valid;
        // Read latitudes and longitudes:
        v_lon.getVar(lon.data());
        v_lat.getVar(lat.data());
        for (int t=0; t<nt; ++t) {
            // Read time record t:
            v.getVar({(size_t)t,0,0}, {1,ny,nx}, xy_slice.data());
            for (int i=0; i<nxy; ++i) {
                if (!std::isnan(xy_slice[i])) {
                    xy_slice_valid.push_back(xy_slice[i]);
                    lon_valid.push_back(lon[i]);
                    lat_valid.push_back(lat[i]);
                }
            }
            // Set up interpolant:
            EarthInterpolant F(lat_valid, lon_valid, xy_slice_valid, (size_t)8, 2.0);
            // Replace NaNs by extrapolated values:
            for (int i=0; i<nxy; ++i) {
                if (std::isnan(xy_slice[i])) {
                    xy_slice[i] = F(lat[i], lon[i]);
                }
            }
            v.putVar({(size_t)t,0,0}, {1,ny,nx}, xy_slice.data());
            // Remove contents of nonan vectors:
            xy_slice_valid.clear();
            lon_valid.clear();
            lat_valid.clear();
        }
    }
    std::cout << "------- END of ROMS_forcing::era_extrapolate_land -------\n";
    std::cout.flush();
}

void ROMS_forcing::extrapolate_land(std::string ROMS_file, std::string ROMS_gridfile,
        std::vector<std::string> variables)
{
    std::println("\n------- ROMS_forcing::extrapolate_land -------");
    std::println("ROMS file:  {}", ROMS_file);
    NcFile nc, nc_grd;
    NcVar v, v_lon, v_lat;
    bool var_on_ugrid, var_on_vgrid;
    size_t nt, ny, nx, nxy;
    //double *lat, *lon, *xy_slice;
    // Open Netcdf files:
    nc_grd.open(ROMS_gridfile, NcFile::read);
    nc.open(ROMS_file, NcFile::FileMode::write);
    nt = nc.getDim("time").getSize();
    for (const auto &var : variables) {
        std::cout << "var = " << var << '\n';
        std::cout.flush();
        v = nc.getVar(var);
        var_on_ugrid = false;
        var_on_vgrid = false;
        // Figure out on which grid this variable lives:
        for (auto &dim : v.getDims()) {
            if (dim.getName() == "xi_u") var_on_ugrid = true;
            if (dim.getName() == "eta_v") var_on_vgrid = true;
        }
        if (var_on_ugrid) {
            // sustr:
            v_lon = nc_grd.getVar("lon_u");
            v_lat = nc_grd.getVar("lat_u");
            ny = nc.getDim("eta_rho").getSize();
            nx = nc.getDim("xi_u").getSize();
        }
        else if (var_on_vgrid) {
            // svstr:
            v_lon = nc_grd.getVar("lon_v");
            v_lat = nc_grd.getVar("lat_v");
            ny = nc.getDim("eta_v").getSize();
            nx = nc.getDim("xi_rho").getSize();
        }
        else {
            // rho grid:
            v_lon = nc_grd.getVar("lon_rho");
            v_lat = nc_grd.getVar("lat_rho");
            ny = nc.getDim("eta_rho").getSize();
            nx = nc.getDim("xi_rho").getSize();
        }
        // Allocate buffers for data, latitudes and longitudes:
        nxy = nx*ny;
        std::vector<double> xy_slice(nxy);
        std::vector<double> lon(nxy), lat(nxy);
        std::vector<double> xy_slice_valid;
        std::vector<double> lon_valid, lat_valid;
        bool has_fillvalue = false;
        double fill_value;
        // Check if this variable has a fill value attribute:
        auto v_attrib = v.getAtts();
        if (v_attrib.contains("_FillValue")) {
            has_fillvalue = true;
            v.getAtt("_FillValue").getValues(&fill_value);
        }
        // Read latitudes and longitudes:
        v_lon.getVar(lon.data());
        v_lat.getVar(lat.data());
        for (int t=0; t<nt; ++t) {
            // Read time record t:
            v.getVar({(size_t)t,0,0}, {1,ny,nx}, xy_slice.data());
            if (has_fillvalue) {
                for (int i=0; i<nxy; ++i) {
                    if (!std::isnan(xy_slice[i]) && xy_slice[i]!=fill_value) {
                        xy_slice_valid.push_back(xy_slice[i]);
                        lon_valid.push_back(lon[i]);
                        lat_valid.push_back(lat[i]);
                    }
                }
            }
            else {
                for (int i=0; i<nxy; ++i) {
                    if (!std::isnan(xy_slice[i])) {
                        xy_slice_valid.push_back(xy_slice[i]);
                        lon_valid.push_back(lon[i]);
                        lat_valid.push_back(lat[i]);
                    }
                }
            }
            // Set up interpolant:
            EarthInterpolant F(lat_valid, lon_valid, xy_slice_valid, (size_t)8, 2.0);
            // Replace NaNs by extrapolated values:
            if (has_fillvalue) {
                for (int i=0; i<nxy; ++i) {
                    if (std::isnan(xy_slice[i]) || xy_slice[i]==fill_value) {
                        xy_slice[i] = F(lat[i], lon[i]);
                    }
                }
            }
            else {
                for (int i=0; i<nxy; ++i) {
                    if (std::isnan(xy_slice[i])) {
                        xy_slice[i] = F(lat[i], lon[i]);
                    }
                }
            }
            v.putVar({(size_t)t,0,0}, {1,ny,nx}, xy_slice.data());
            // Remove contents of nonan vectors:
            xy_slice_valid.clear();
            lon_valid.clear();
            lat_valid.clear();
        }
    }
    std::cout << "------- END of ROMS_forcing::extrapolate_land -------\n";
    std::cout.flush();
}