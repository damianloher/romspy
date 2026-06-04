#include <chrono>
#include <print>
#include "roms_forcing_utils.h"

ROMS_forcing::ROMS_forcing(std::string ROMS_setup, std::string ROMS_gridfile, 
        std::string outdir, bool verbose)
{
    std::println("ROMS_forcing class constructor: ROMS_gridfile = {}", ROMS_gridfile);
    this->ROMS_setup = ROMS_setup;
    this->ROMS_gridfile = ROMS_gridfile;
    this->output_dir = outdir;
    this->verbose = verbose;
    std::cout.flush();
}

// Creates a Netcdf file containing relative humidty based of t2m and d2m ERA5 data. The
// ERA5 data is assumed to be on a ROMS grid.
//
// Arguments:
// - outfile: output file (will be overwritten if it exists already)
// - year: year of the data, used only in the name of the output file
// - t2m_file: Netcdf file containing temperature at 2m (Celsius) on ROMS grid
// - d2m_file: Netcdf file containing dew point temperature at 2m (K) on ROMS grid
void ROMS_forcing::rel_humidity(std::string outfile, int year, std::string t2m_file, std::string d2m_file) {
    std::cout << "\n------- ROMS_forcing::Relative humidity -------";
    std::cout << "year: " << year << '\n';
    NcFile nc_t2m, nc_d2m, nc_rhum;
    NcVar v_d2m, v_t2m, v_qair, v_lon, v_lat;
    std::map<std::string,int> dim_size;
    std::multimap<std::string,NcDim> ncdims;
    int ndims;
    // Open t2m and d2m files:
    nc_t2m.open(t2m_file, NcFile::read);
    nc_d2m.open(d2m_file, NcFile::read);
    // Store sizes of all dimensions in map "dim_size":
    ncdims = nc_t2m.getDims();
    ndims = nc_t2m.getDimCount();
    for (auto it = ncdims.begin(); it != ncdims.end(); it++) {
        dim_size[it->second.getName()] = it->second.getSize();
    }
    // Size of xy slice:
    size_t nxy = dim_size["eta_rho"]*dim_size["xi_rho"];

    // Create output file: if the file exists already, it is overwritten
    std::cout << "output file: " << outfile << '\n';
    NcFile nc_out(outfile, NcFile::FileMode::replace,
                NcFile::FileFormat::nc4);
    // Create dimensions:
    for (const auto& [dim, dsize] : dim_size) {
        nc_out.addDim(dim, dsize);
    }
    const char *attval = new char[128];
    char *attval_buf = new char[128];
    // Create variable "time":
    NcVar v_time  = nc_out.addVar("time", "double", "time");
    attval = "time";
    v_time.putAtt("standard_name", NcType::nc_CHAR, strlen(attval), attval);
    attval = "proleptic_gregorian";
    v_time.putAtt("calendar", NcType::nc_CHAR, strlen(attval), attval);
    attval = "T";
    v_time.putAtt("axis", NcType::nc_CHAR, strlen(attval), attval);
    sprintf(attval_buf, "hours since %i-01-01 00:00:00", year);
    v_time.putAtt("units", NcType::nc_CHAR, strlen(attval_buf), attval_buf);
    // Copy times from t2m file:
    v_t2m = nc_t2m.getVar("time");
    double *time_buf = new double[dim_size["time"]];
    v_t2m.getVar(time_buf);
    v_time.putVar(time_buf);
    delete[] time_buf;
    // Longitudes of rho points:
    std::vector<std::string> dims_2d = {"eta_rho", "xi_rho"};
    v_lon  = nc_out.addVar("lon_rho", "double", dims_2d);
    attval = "longitude";
    v_lon.putAtt("standard_name", NcType::nc_CHAR, strlen(attval), attval);
    attval = "longitude of RHO-points";
    v_lon.putAtt("long_name", NcType::nc_CHAR, strlen(attval), attval);
    attval = "degree_east";
    v_lon.putAtt("units", NcType::nc_CHAR, strlen(attval), attval);
    attval = "Lon";
    v_lon.putAtt("_CoordinateAxisType", NcType::nc_CHAR, strlen(attval), attval);
    // Allocate buffer for one time slice of data:
    double *dat_buf = new double[nxy];
    // Copy longitudes from t2m_file:
    v_t2m = nc_t2m.getVar("lon_rho");
    v_t2m.getVar(dat_buf);
    v_lon.putVar(dat_buf);
    // Latitudes of rho points:
    v_lat  = nc_out.addVar("lat_rho", "double", dims_2d);
    attval = "latitude";
    v_lat.putAtt("standard_name", NcType::nc_CHAR, strlen(attval), attval);
    attval = "latitude of RHO-points";
    v_lat.putAtt("long_name", NcType::nc_CHAR, strlen(attval), attval);
    attval = "degree_north";
    v_lat.putAtt("units", NcType::nc_CHAR, strlen(attval), attval);
    attval = "Lat";
    v_lat.putAtt("_CoordinateAxisType", NcType::nc_CHAR, strlen(attval), attval);
    // Copy latitudes from t2m_file:
    v_t2m = nc_t2m.getVar("lat_rho");
    v_t2m.getVar(dat_buf);
    v_lat.putVar(dat_buf);
    // Create variable "Qair":
    v_qair  = nc_out.addVar("Qair", "float", {"time", "eta_rho", "xi_rho"});
    attval = "relative_humidity";
    v_qair.putAtt("standard_name", NcType::nc_CHAR, strlen(attval), attval);
    attval = "Near-Surface Relative Humidity";
    v_qair.putAtt("long_name", NcType::nc_CHAR, strlen(attval), attval);
    attval = "fraction";
    v_qair.putAtt("units", NcType::nc_CHAR, strlen(attval), attval);
    attval = "lon_rho lat_rho";
    v_qair.putAtt("coordinates", NcType::nc_CHAR, strlen(attval), attval);

    // Compute and write relative humidity:
    // References:
    // [1] ECMWF IFS Documentation Cy41r2, part IV Physical Processes
    // [2] https://confluence.ecmwf.int/display/CKB/ERA5%3A+data+documentation
    //     (guidelines, item 10)
    // Constants used in [1, eq. 7.5] for saturation over water, see [1, p. 95]:
    //double a1 = 611.21;  // Pa
    double a3 = 17.502;  // has no unit
    double a4 = 32.19;   // K
    double T0 = 273.16;  // K
    std::vector<size_t> startp, countp;
    startp.push_back(0);
    startp.push_back(0);
    startp.push_back(0);
    countp.push_back(1);
    countp.push_back(dim_size["eta_rho"]);
    countp.push_back(dim_size["xi_rho"]);
    v_t2m = nc_t2m.getVar("t2m");
    v_d2m = nc_d2m.getVar("d2m");
    double *dat_buf2 = new double[nxy];
    double *rel_hum = new double[nxy];
    double esat, esat_d;
    for (size_t t=0; t<dim_size["time"]; t++) {
        startp[0] = t;
        v_t2m.getVar(startp,countp,dat_buf);  // t2m, Celsius
        v_d2m.getVar(startp,countp,dat_buf2);  // d2m, K
        for (size_t i=0; i<nxy; i++) {
            // Saturation vapor pressures [1, eq. 7.5]:
            esat = std::exp(a3*dat_buf[i]/(dat_buf[i]+T0-a4));
            esat_d = std::exp(a3*(dat_buf2[i]-T0)/(dat_buf2[i]-a4));
            // Relative humidity as a fraction [2]:
            rel_hum[i] = esat_d/esat;
        }
        // Write relative humidity to output file:
        v_qair.putVar(startp,countp,rel_hum);
    }
    // Delete allocated buffers:
    delete[] dat_buf, dat_buf2, rel_hum;
    std::cout << "------- END of ROMS_forcing::Relative humidity -------\n";
    std::cout.flush();
}

void ROMS_forcing::wind_stress(std::string outfile, int year,
        std::string era5_ewss_file, std::string era5_nsss_file,
        std::string ref_file, std::string time_dim_ref, int tstep_day)
{
    std::cout << "--- Wind stress ---\n";
    std::cout << "year: " << year << '\n';
    NcFile nc_grd, nc_ewss, nc_nsss, nc_ref;
    NcVar v_angle, v_sustr, v_svstr, v_ew, v_ns, v_lon, v_lat,
        v_ref, v_grid;
    std::map<std::string,int> dim_size;
    std::multimap<std::string,NcDim> ncdims;
    int ndims;
    // Open ROMS grid file:
    nc_grd.open(ROMS_gridfile, NcFile::read);
    // Store sizes of dimensions in ROMS grid file, in map "dim_size":
    ncdims = nc_grd.getDims();
    ndims = nc_grd.getDimCount();
    for (auto it = ncdims.begin(); it != ncdims.end(); it++) {
        dim_size[it->second.getName()] = it->second.getSize();
    }
    // Sizes of xy slices on various grids:
    size_t nxy = dim_size["eta_rho"]*dim_size["xi_rho"];
    size_t nxy_u = dim_size["eta_u"]*dim_size["xi_u"];
    size_t nxy_v = dim_size["eta_v"]*dim_size["xi_v"];
    // Allocate buffers:
    double *dat_buf1 = new double[nxy];
    double *dat_buf2 = new double[nxy];
    double *dat_buf_u = new double[nxy_u];
    double *dat_buf_v = new double[nxy_v];
    double *angle = new double[nxy];
    double *cos_angle = new double[nxy];
    double *sin_angle = new double[nxy];
    // Read angles and compute sin and cos of them:
    v_angle = nc_grd.getVar("angle");
    v_angle.getVar(angle);
    for (int i=0; i<nxy; i++) {
        cos_angle[i] = std::cos(angle[i]);
        sin_angle[i] = std::sin(angle[i]);
    }
    // Open reference file (used for reading times):
    nc_ref.open(ref_file, NcFile::read);

    // Create output file: if the file exists already, it is overwritten
    std::cout << "output file: " << outfile << '\n';
    std::cout << std::flush;
    NcFile nc_out(outfile, NcFile::FileMode::replace,
                NcFile::FileFormat::nc4);
    // Create dimensions present in the ROMS grid file:
    for (const auto& [dim, dsize] : dim_size) {
        nc_out.addDim(dim, dsize);
    }
    // Create time dimension:
    NcDim d_time = nc_ref.getDim(time_dim_ref);
    dim_size["sms_time"] = d_time.getSize();
    nc_out.addDim("sms_time", dim_size["sms_time"]);
    const char *attval = new char[128];
    char *attval_buf = new char[128];
    double attval_d;
    // Create variable "sms_time":
    v_ref = nc_ref.getVar(time_dim_ref);
    NcVar v_time  = nc_out.addVar("sms_time", "double", "sms_time");
    attval = "sms_time";
    v_time.putAtt("long_name", NcType::nc_CHAR, strlen(attval), attval);
    attval_d = 365.0;
    if (std::chrono::year(year).is_leap()) {
        attval_d = 366.0;
    }
    v_time.putAtt("cycle_length", NcType::nc_DOUBLE, attval_d);
    attval = "proleptic_gregorian";
    v_time.putAtt("calendar", NcType::nc_CHAR, strlen(attval), attval);
    attval = "T";
    v_time.putAtt("axis", NcType::nc_CHAR, strlen(attval), attval);
    // Copy "units" attribute from ref file:
    v_ref.getAtt("units").getValues(attval_buf);
    v_time.putAtt("units", NcType::nc_CHAR, strlen(attval_buf), attval_buf);
    // Copy times from ref file:
    double *time_buf = new double[dim_size["time"]];
    v_ref.getVar(time_buf);
    v_time.putVar(time_buf);
    delete[] time_buf;
    // Longitudes of rho points:
    std::vector<std::string> dims_2d = {"eta_rho", "xi_rho"};
    v_lon  = nc_out.addVar("lon_rho", "double", dims_2d);
    attval = "longitude";
    v_lon.putAtt("standard_name", NcType::nc_CHAR, strlen(attval), attval);
    attval = "longitude of RHO-points";
    v_lon.putAtt("long_name", NcType::nc_CHAR, strlen(attval), attval);
    attval = "degree_east";
    v_lon.putAtt("units", NcType::nc_CHAR, strlen(attval), attval);
    attval = "Lon";
    v_lon.putAtt("_CoordinateAxisType", NcType::nc_CHAR, strlen(attval), attval);
    // Copy longitudes from ROMS grid file:
    v_grid = nc_grd.getVar("lon_rho");
    v_grid.getVar(dat_buf1);
    v_lon.putVar(dat_buf1);
    // Latitudes of rho points:
    v_lat  = nc_out.addVar("lat_rho", "double", dims_2d);
    attval = "latitude";
    v_lat.putAtt("standard_name", NcType::nc_CHAR, strlen(attval), attval);
    attval = "latitude of RHO-points";
    v_lat.putAtt("long_name", NcType::nc_CHAR, strlen(attval), attval);
    attval = "degree_north";
    v_lat.putAtt("units", NcType::nc_CHAR, strlen(attval), attval);
    attval = "Lat";
    v_lat.putAtt("_CoordinateAxisType", NcType::nc_CHAR, strlen(attval), attval);
    // Copy latitudes from ROMS grid file:
    v_grid = nc_grd.getVar("lat_rho");
    v_grid.getVar(dat_buf1);
    v_lat.putVar(dat_buf1);
    // Create variable "sustr":
    v_sustr  = nc_out.addVar("sustr", "float", {"sms_time", "eta_u", "xi_u"});
    attval = "surface u-momentum stress";
    v_sustr.putAtt("long_name", NcType::nc_CHAR, strlen(attval), attval);
    attval = "Newton meter-2";
    v_sustr.putAtt("units", NcType::nc_CHAR, strlen(attval), attval);
    attval = "ERA5 ewss and nsss data";
    v_sustr.putAtt("source", NcType::nc_CHAR, strlen(attval), attval);
    // Create variable "svstr":
    v_svstr  = nc_out.addVar("svstr", "float", {"sms_time", "eta_v", "xi_v"});
    attval = "surface v-momentum stress";
    v_svstr.putAtt("long_name", NcType::nc_CHAR, strlen(attval), attval);
    attval = "Newton meter-2";
    v_svstr.putAtt("units", NcType::nc_CHAR, strlen(attval), attval);
    attval = "ERA5 ewss and nsss data";
    v_svstr.putAtt("source", NcType::nc_CHAR, strlen(attval), attval);

    // Scaling factor for ERA5 ewss and nsss data:
    double stimei = tstep_day/(24*3600);  // s-1
    // Initialize start and count vectors for reading/writing data:
    std::vector<size_t> startp, countp_rho, countp_u, countp_v;
    startp.push_back(0);
    startp.push_back(0);
    startp.push_back(0);
    countp_rho.push_back(1);
    countp_rho.push_back(dim_size["eta_rho"]);
    countp_rho.push_back(dim_size["xi_rho"]);
    countp_u.push_back(1);
    countp_u.push_back(dim_size["eta_u"]);
    countp_u.push_back(dim_size["xi_u"]);
    countp_v.push_back(1);
    countp_v.push_back(dim_size["eta_v"]);
    countp_v.push_back(dim_size["xi_v"]);
    nc_ewss.open(era5_ewss_file, NcFile::read);
    nc_nsss.open(era5_nsss_file, NcFile::read);
    v_ew = nc_ewss.getVar("ewss");
    v_ns = nc_nsss.getVar("nsss");
    double *ew_rot = new double[nxy];
    double *ns_rot = new double[nxy];
    size_t nx_u = dim_size["xi_u"];
    size_t nx_v = dim_size["xi_v"];
    for (size_t t=0; t<dim_size["sms_time"]; t++) {
        startp[0] = t;
        v_ew.getVar(startp,countp_rho,dat_buf1);  // Newton m-2 s
        v_ns.getVar(startp,countp_rho,dat_buf2);  // Newton m-2 s
        // Rotate to ROMS rho grid and convert to N m-2:
        for (size_t i=0; i<nxy; i++) {
            ew_rot[i] = stimei * (cos_angle[i]*dat_buf1[i] + sin_angle[i]*dat_buf2[i]);
            ns_rot[i] = stimei * (cos_angle[i]*dat_buf2[i] - sin_angle[i]*dat_buf1[i]);
        }
        // Interpolate sustr to u grid:
        for (size_t i=0; i<dim_size["eta_u"]; i++) {
            for (size_t j=0; j<nx_u; j++) {
                dat_buf_u[i*nx_u+j] = 0.5*(ew_rot[i*nx_u+j]+ew_rot[i*nx_u+j+1]);
            }
        }
        // Interpolate svstr to v grid:
        for (size_t i=0; i<dim_size["eta_v"]; i++) {
            for (size_t j=0; j<nx_v; j++) {
                dat_buf_v[i*nx_v+j] = 0.5*(ns_rot[i*nx_v+j]+ns_rot[(i+1)*nx_v+j]);
            }
        }
        // Write data to output file:
        v_sustr.putVar(startp,countp_u,dat_buf_u);
        v_svstr.putVar(startp,countp_v,dat_buf_v);
    }
}
