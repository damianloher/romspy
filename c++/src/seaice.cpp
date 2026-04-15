#include <chrono>
#include <limits>
#include <netcdf>
#include <print>
#include <stdexcept>
#include <vector>
// Enable range checking of boost multiarrays by defining the following:
#define BOOST_DISABLE_ASSERTS
#include <boost/multi_array.hpp>
#include "roms_forcing_utils.h"


void ROMS_forcing::make_seaice_correction(std::vector<std::string> frc_files, int year, std::string seaice_file,
    std::string snowice_file, int time_res_frc_h, bool verbose)
{
    std::println("\n------- START of ROMS_forcing::make_seaice_correction -------");
    std::println("daily sea-ice file:  {}", seaice_file);
    std::println("daily snow-ice file: {}", snowice_file);
    std::cout.flush();
    NcFile nc_seaice, nc_snowice;
    NcVar var;
    int nr_days_per_year;
    float cycle_length;
    // Number of forcing time records per day:
    float day_per_trec = time_res_frc_h / 24.0f;
    std::vector<float> seaice_time_frc, seaice_time_dat;
    if (std::chrono::year(year).is_leap()) {
        nr_days_per_year= 366;
        cycle_length = 366.0;
    } else {
        nr_days_per_year= 365;
        cycle_length = 365.0;
    }
    // Set time values (in days) of the forcing data:
    // these values are used for variable ice_time and for
    // interpolating the daily sea ice data to the time resolution
    // of the ROMS forcing
    for (float f=0.5*day_per_trec; f<nr_days_per_year+0.1*day_per_trec; f+=day_per_trec) {
        seaice_time_frc.push_back(f);
    }
    // Time values (in days) of the sea and snow ice daily data:
    for (float f=0.5; f<nr_days_per_year; f+=1.0) seaice_time_dat.push_back(f);

    // Sea ice related data:
    std::println("prepare sea ice data");
    nc_seaice.open(seaice_file, NcFile::read);
    // Store sizes of all dimensions in map "dim_size":
    std::map<std::string,int> dim_size;
    std::multimap<std::string,NcDim> ncdims;
    int ndims;
    ncdims = nc_seaice.getDims();
    ndims = nc_seaice.getDimCount();
    for (auto it = ncdims.begin(); it != ncdims.end(); ++it) {
        dim_size[it->second.getName()] = it->second.getSize();
    }
    unsigned long int nt_seaice = dim_size["ice_time"];
    unsigned long int nx = dim_size["xi_rho"];
    unsigned long int ny = dim_size["eta_rho"];
    unsigned long int nx_sustr, ny_svstr;
    // Daily sea ice concentration for 366d:
    typedef boost::multi_array<float,2> array_2d;
    typedef boost::multi_array<float,3> array_3d;
    typedef array_3d::index index;
    // Set up sea ice related data:
    array_3d seaice_daily(boost::extents[nr_days_per_year][ny][nx]);
    array_3d shflux_seaice_daily(boost::extents[nr_days_per_year][ny][nx]);
    array_3d swflux_seaice_daily(boost::extents[nr_days_per_year][ny][nx]);
    std::vector<size_t> start, count;
    // Check if the daily sea ice data needs to be extended (i.e. if Feb 29
    // is not present for a leap year):
    if (nt_seaice == 365 && nr_days_per_year== 366) {
        // Extend data to 366 days:
        var = nc_seaice.getVar("seaice");
        // Read January and February (time records 0 to 58):
        start = {0,0,0};
        count = {59,ny,nx};
        var.getVar(start, count, seaice_daily.data());
        // Read time records 59 to 364 into slots 60 to 365 of seaice_daily:
        start = {59,0,0};
        count = {306,ny,nx};
        var.getVar(start, count, &seaice_daily[60][0][0]);
        // Compute data for Feb 29 (slot 59 of seaice_daily):
        for (int j=0; j<ny; ++j)
            for (int i=0; i<nx; ++i)
                seaice_daily[59][j][i] = 0.5*(seaice_daily[58][j][i]+seaice_daily[60][j][i]);
        // Surface net heat flux:
        var = nc_seaice.getVar("shflux");
        // Read January and February (time records 0 to 58):
        start = {0,0,0};
        count = {59,ny,nx};
        var.getVar(start, count, shflux_seaice_daily.data());
        // Read time records 59 to 364 into slots 60 to 365 of seaice_daily:
        start = {59,0,0};
        count = {306,ny,nx};
        var.getVar(start, count, &shflux_seaice_daily[60][0][0]);
        for (int j=0; j<ny; ++j)
            for (int i=0; i<nx; ++i)
                shflux_seaice_daily[59][j][i] = 0.5*(shflux_seaice_daily[58][j][i]
                            + shflux_seaice_daily[60][j][i]);
        // Surface fresh water flux:
        var = nc_seaice.getVar("swflux");
        // Read January and February (time records 0 to 58):
        start = {0,0,0};
        count = {59,ny,nx};
        var.getVar(start, count, swflux_seaice_daily.data());
        // Read time records 59 to 364 into slots 60 to 365 of seaice_daily:
        start = {59,0,0};
        count = {306,ny,nx};
        var.getVar(start, count, &swflux_seaice_daily[60][0][0]);
        for (int j=0; j<ny; ++j)
            for (int i=0; i<nx; ++i)
                swflux_seaice_daily[59][j][i] = 0.5*(swflux_seaice_daily[58][j][i]
                            + swflux_seaice_daily[60][j][i]);
    } else if (nt_seaice == nr_days_per_year) {
        // Data in Netcdf file already has the correct number of time records:
        var = nc_seaice.getVar("seaice");
        var.getVar(seaice_daily.data());
        var = nc_seaice.getVar("shflux");
        var.getVar(shflux_seaice_daily.data());
        var = nc_seaice.getVar("swflux");
        var.getVar(swflux_seaice_daily.data());
    } else if (nt_seaice == 366 && nr_days_per_year== 365) {
        auto msg = std::format("not implemented: {} sea ice time records found for non-leap year", nt_seaice);
        throw std::runtime_error(msg);
    } else {
        auto msg = std::format("number of time records of sea ice file must be 365 or 366 (found {})", nt_seaice);
        throw std::runtime_error(msg);
    }
    nc_seaice.close();

    // Set up snow ice related data:
    std::println("prepare snow ice data");
    std::cout.flush();
    nc_snowice.open(snowice_file, NcFile::read);
    // Store sizes of all dimensions in map "dim_size_snow":
    std::map<std::string,int> dim_size_snow;
    ncdims = nc_snowice.getDims();
    ndims = nc_snowice.getDimCount();
    for (auto it = ncdims.begin(); it != ncdims.end(); ++it) {
        dim_size_snow[it->second.getName()] = it->second.getSize();
    }
    unsigned long int nt_snowice = dim_size_snow["ice_time"];
    unsigned long int nx_snowice = dim_size_snow["xi_rho"];
    unsigned long int ny_snowice = dim_size_snow["eta_rho"];
    if (nx != nx_snowice) {
        auto msg = std::format("size of dimension 'xi_rho' differs between seaice and snowice file: {} vs {}", nx, nx_snowice);
        throw std::runtime_error(msg);
    }
    if (ny != ny_snowice) {
        auto msg = std::format("size of dimension 'eta_rho' differs between seaice and snowice file: {} vs {}", ny, ny_snowice);
        throw std::runtime_error(msg);
    }
    array_3d snowice_daily(boost::extents[nr_days_per_year][ny][nx]);
    array_3d shflux_snowice_daily(boost::extents[nr_days_per_year][ny][nx]);
    array_3d swflux_snowice_daily(boost::extents[nr_days_per_year][ny][nx]);
    // Check if the daily snow ice data needs to be extended (i.e. if Feb 29
    // is not present for a leap year):
    if (nt_snowice == 365 && nr_days_per_year== 366) {
        // Extend data to 366 days:
        var = nc_snowice.getVar("seaice");
        // Read January and February (time records 0 to 58):
        start = {0,0,0};
        count = {59,ny,nx};
        var.getVar(start, count, snowice_daily.data());
        // Read time records 59 to 364 into slots 60 to 365 of snowice_daily:
        start = {59,0,0};
        count = {306,ny,nx};
        var.getVar(start, count, &snowice_daily[60][0][0]);
        // Compute data for Feb 29 (slot 59 of snowice_daily):
        for (int j=0; j<ny; ++j)
            for (int i=0; i<nx; ++i)
                snowice_daily[59][j][i] = 0.5*(snowice_daily[58][j][i]+snowice_daily[60][j][i]);
        // Surface net heat flux:
        var = nc_snowice.getVar("shflux");
        // Read January and February (time records 0 to 58):
        start = {0,0,0};
        count = {59,ny,nx};
        var.getVar(start, count, shflux_snowice_daily.data());
        // Read time records 59 to 364 into slots 60 to 365 of snowice_daily:
        start = {59,0,0};
        count = {306,ny,nx};
        var.getVar(start, count, &shflux_snowice_daily[60][0][0]);
        for (int j=0; j<ny; ++j)
            for (int i=0; i<nx; ++i)
                shflux_snowice_daily[59][j][i] = 0.5*(shflux_snowice_daily[58][j][i]
                            + shflux_snowice_daily[60][j][i]);
        // Surface fresh water flux:
        var = nc_snowice.getVar("swflux");
        // Read January and February (time records 0 to 58):
        start = {0,0,0};
        count = {59,ny,nx};
        var.getVar(start, count, swflux_snowice_daily.data());
        // Read time records 59 to 364 into slots 60 to 365 of snowice_daily:
        start = {59,0,0};
        count = {306,ny,nx};
        var.getVar(start, count, &swflux_snowice_daily[60][0][0]);
        for (int j=0; j<ny; ++j)
            for (int i=0; i<nx; ++i)
                swflux_snowice_daily[59][j][i] = 0.5*(swflux_snowice_daily[58][j][i]
                            + swflux_snowice_daily[60][j][i]);
    } else if (nt_snowice == nr_days_per_year) {
        // Data in files already has the correct number of time records:
        var = nc_snowice.getVar("seaice");
        var.getVar(snowice_daily.data());
        var = nc_snowice.getVar("shflux");
        var.getVar(shflux_snowice_daily.data());
        var = nc_snowice.getVar("swflux");
        var.getVar(swflux_snowice_daily.data());
    } else if (nt_snowice == 366 && nr_days_per_year== 365) {
        auto msg = std::format("not implemented: {} snow ice time records found for non-leap year", nt_snowice);
        throw std::runtime_error(msg);
    } else {
        auto msg = std::format("number of time records of snow ice file must be 365 or 366 (found {})", nt_snowice);
        throw std::runtime_error(msg);
    }
    nc_snowice.close();

    array_2d seaice_curr(boost::extents[ny][nx]);
    array_2d snowice_curr(boost::extents[ny][nx]);
    array_2d shflux_seaice_curr(boost::extents[ny][nx]);
    array_2d swflux_seaice_curr(boost::extents[ny][nx]);
    array_2d shflux_snowice_curr(boost::extents[ny][nx]);
    array_2d swflux_snowice_curr(boost::extents[ny][nx]);
    array_2d albedo(boost::extents[ny][nx]);
    array_2d swrad_frc(boost::extents[ny][nx]);
    array_2d swflux_frc(boost::extents[ny][nx]);
    array_2d shflux_frc(boost::extents[ny][nx]);
    array_2d shflux_ice_frc(boost::extents[ny][nx]);
    array_2d swrad_ice_frc(boost::extents[ny][nx]);
    int i1, i2, t1, t2, idx;
    float cff, time_factor, tmp, tmp1, tmp2;
    float shflux_ice_lat, shflux_ice_melt, shflux_melt, shflux_ice_freeze, shflux_ice_sen;
    size_t nt_frc;
    bool sustr_on_ugrid, svstr_on_vgrid;
    std::multimap<std::string, NcVar> vars_frc;
    std::vector<NcDim> dim_vec;
    NcFile nc_frc;
    NcVar vobj, v_seaice, v_swrad_frc, v_shflux_frc, v_swflux_frc, v_sustr, v_svstr;
    std::string time_unit;
    std::string del = " ";
    std::string::size_type pos1, pos2, pos3;
    int tunit_year, tunit_month, tunit_day;
    std::chrono::year_month_day start_data_time, start_frc_time;
    const char *attval = new char[128];
    char *attval_buf = new char[128];
    t1 = 0;
    t2 = 0;
    // Loop for computing nansum(shflux_melt) and nansum(shflux_ice):
    double nansum_shflux_melt = 0.0;
    double nansum_shflux_ice = 0.0;
    std::println("compute nansums");
    std::cout.flush();
    for (auto f: frc_files) {
        nc_frc.open(f, NcFile::write);
        nt_frc = nc_frc.getDim("time").getSize();
        //std::println("nt_frc = {}", nt_frc);
        t1 = t2;
        t2 += nt_frc;
        //std::println("t1 = {}", t1);
        //std::println("t2 = {}", t2);
        // Loop over forcing times (in seaice_time_frc) and interpolate daily data 
        // variables to the forcing times:
        v_shflux_frc = nc_frc.getVar("shflux");
        for (unsigned long t=0; t<nt_frc; ++t) {
            idx = 0;
            while (idx<nr_days_per_year && seaice_time_dat[idx]<=seaice_time_frc[t1+t]) ++idx;
            if (idx == 0) {
                // Before first daily time record:
                i1 = 0;
                i2 = nr_days_per_year-1;
                cff = 0.5 - seaice_time_frc[t1+t];
            } else if (idx == nr_days_per_year) {
                // At or past last daily time record:
                i1 = nr_days_per_year-1;
                i2 = 0;
                cff = seaice_time_frc[t1+t] - nr_days_per_year-0.5;
            } else {
                // seaice_time_dat[idx-1] <= seaice_time_frc[t1+t] < seaice_time_dat[idx]:
                i1 = idx-1;
                i2 = idx;
                cff = seaice_time_frc[t1+t] - seaice_time_dat[i1];
            }
            // Interpolate daily variables to current forcing time:
            for (int j=0; j<ny; ++j)
                for (int i=0; i<nx; ++i) {
                    seaice_curr[j][i] = seaice_daily[i1][j][i] + cff*(seaice_daily[i2][j][i]-seaice_daily[i1][j][i]);
                    shflux_seaice_curr[j][i] = shflux_seaice_daily[i1][j][i] + cff*(
                        shflux_seaice_daily[i2][j][i]-shflux_seaice_daily[i1][j][i]);
                    shflux_snowice_curr[j][i] = shflux_snowice_daily[i1][j][i] + cff*(
                        shflux_snowice_daily[i2][j][i]-shflux_snowice_daily[i1][j][i]);
                }
            // Heat flux:
            v_shflux_frc.getVar({t,0,0}, {1,ny,nx}, shflux_frc.data());
            for (int j=0; j<ny; ++j)
                for (int i=0; i<nx; ++i) {
                    tmp = shflux_frc[j][i];  // shflux_atm in Matlab code
                    shflux_ice_frc[j][i] = tmp*seaice_curr[j][i];  // heat flux into sea ice
                    shflux_frc[j][i] = tmp*(1.0-seaice_curr[j][i]);  // heat flux that goes directly into the ocean
                }
            for (int j=0; j<ny; ++j)
                for (int i=0; i<nx; ++i) {
                    shflux_ice_lat = 1.3*shflux_seaice_curr[j][i]; // latent heat flux associated with melting and freezing
                    shflux_ice_melt = shflux_ice_lat>0 ? 0.0f : shflux_ice_lat;
                    //shflux_ice_melt(shflux_ice_lat>0)=0;
                    //shflux_ice_freeze=shflux_ice_lat;
                    //shflux_ice_freeze(shflux_ice_lat<0)=0;
                    shflux_ice_freeze = shflux_ice_lat<0 ? 0.0f : shflux_ice_lat;
                    // Ice melts from below and freezes from above:
                    // i.e. we add the latent heat of freezing to the total ice
                    // flux and also the short wave flux on top of the ice. 
                    // The net is sensible heat into the ice the acts to
                    // cool or warm the ice layer (ice heat capacity). This heat
                    // stored in the ice is released upon melting to the ocean (will
                    // cool the ocean if negative).
                    shflux_ice_frc[j][i] += shflux_ice_freeze+swrad_ice_frc[j][i];
                    shflux_melt = shflux_ice_melt + shflux_snowice_curr[j][i];  // total latent heat from sea ice and snow melt
                    if (!std::isnan(shflux_ice_frc[j][i])) nansum_shflux_ice += shflux_ice_frc[j][i];
                    if (!std::isnan(shflux_melt)) nansum_shflux_melt += shflux_melt;
                }
        }
        nc_frc.close();
    }

    // Loop over files and correct fluxes:
    t1 = 0;
    t2 = 0;
    float albedo_ice = 0.85;
    float albedo_ocean = 0.06;
    float seaice_min = std::numeric_limits<float>::max();
    float seaice_max = std::numeric_limits<float>::min();
    for (auto f: frc_files) {
        std::println("process frc file: {}", f);
        std::cout.flush();
        nc_frc.open(f, NcFile::write);
        nt_frc = nc_frc.getDim("time").getSize();
        // Add dimension "ice_time" and 2 variables, if they are not present:
        t1 = t2;
        t2 += nt_frc;
        if (verbose) {
            std::println("t1 = {}", t1);
            std::println("t2 = {}", t2);
        }
        std::cout.flush();
        // Create dimension "ice_time" if it is not present already:
        if (!nc_frc.getDims().contains("ice_time")) {
            nc_frc.addDim("ice_time", nt_frc);
        }
        // Create variable "ice_time" if not present already:
        vars_frc = nc_frc.getVars();
        if (!vars_frc.contains("ice_time")) {
            //std::println("add variable \"ice_time\"");
            vobj = nc_frc.addVar("ice_time", "float", "ice_time");
            attval = "ice_time";
            vobj.putAtt("long_name", NcType::nc_CHAR, strlen(attval), attval);
            sprintf(attval_buf, "days since %i-01-01 00:00", year);
            vobj.putAtt("units", NcType::nc_CHAR, strlen(attval_buf), attval_buf);
            vobj.putAtt("cycle_length", NcType::nc_FLOAT, cycle_length);
            attval = "means";
            vobj.putAtt("climatological", NcType::nc_CHAR, strlen(attval), attval);
            vobj.putVar({0}, {nt_frc}, &seaice_time_frc[t1]);
        }
        // Create variable "seaice" if not present already:
        if (!vars_frc.contains("seaice")) {
            //std::println("add variable \"seaice\"");
            v_seaice = nc_frc.addVar("seaice", "float", {"ice_time","y","x"});
            v_seaice.setFill(true, -1.0e3f);
            attval = "Seaice fraction";
            v_seaice.putAtt("long_name", NcType::nc_CHAR, strlen(attval), attval);
            attval = "-";
            v_seaice.putAtt("units", NcType::nc_CHAR, strlen(attval), attval);
            attval = seaice_file.c_str();
            v_seaice.putAtt("seaice_corr", NcType::nc_CHAR, strlen(attval), attval);
        } else {
            v_seaice = nc_frc.getVar("seaice");
        }
        // Loop over forcing times (in seaice_time_frc) and interpolate daily data 
        // variables to the forcing times:
        v_swrad_frc = nc_frc.getVar("swrad");
        v_shflux_frc = nc_frc.getVar("shflux");
        v_swflux_frc = nc_frc.getVar("swflux");
        //std::println("apply seaice correction");
        for (unsigned long t=0; t<nt_frc; ++t) {
            idx = 0;
            while (idx<nr_days_per_year && seaice_time_dat[idx]<=seaice_time_frc[t1+t]) ++idx;
            if (idx == 0) {
                // Before first daily time record:
                i1 = 0;
                i2 = nr_days_per_year-1;
                cff = 0.5 - seaice_time_frc[t1+t];
            } else if (idx == nr_days_per_year) {
                // At or past last daily time record:
                i1 = nr_days_per_year-1;
                i2 = 0;
                cff = seaice_time_frc[t1+t] - nr_days_per_year-0.5;
            } else {
                // seaice_time_dat[idx-1] <= seaice_time_frc[t1+t] < seaice_time_dat[idx]:
                i1 = idx-1;
                i2 = idx;
                cff = seaice_time_frc[t1+t] - seaice_time_dat[i1];
            }
            // Interpolate daily variables to current forcing time:
            seaice_min = std::numeric_limits<float>::max();
            seaice_max = std::numeric_limits<float>::min();
            for (int j=0; j<ny; ++j)
                for (int i=0; i<nx; ++i) {
                    seaice_curr[j][i] = seaice_daily[i1][j][i] + cff*(seaice_daily[i2][j][i]-seaice_daily[i1][j][i]);
                    //if (seaice_curr[j][i] < seaice_min) seaice_min = seaice_curr[j][i];
                    //if (seaice_curr[j][i] > seaice_max) seaice_max = seaice_curr[j][i];
                    snowice_curr[j][i] = snowice_daily[i1][j][i] + cff*(snowice_daily[i2][j][i]-snowice_daily[i1][j][i]);
                    shflux_seaice_curr[j][i] = shflux_seaice_daily[i1][j][i] + cff*(
                        shflux_seaice_daily[i2][j][i]-shflux_seaice_daily[i1][j][i]);
                    swflux_seaice_curr[j][i] = swflux_seaice_daily[i1][j][i] + cff*(
                        swflux_seaice_daily[i2][j][i]-swflux_seaice_daily[i1][j][i]);
                    shflux_snowice_curr[j][i] = shflux_snowice_daily[i1][j][i] + cff*(
                        shflux_snowice_daily[i2][j][i]-shflux_snowice_daily[i1][j][i]);
                    swflux_snowice_curr[j][i] = swflux_snowice_daily[i1][j][i] + cff*(
                        swflux_snowice_daily[i2][j][i]-swflux_snowice_daily[i1][j][i]);
                }
            // Write current sea ice fraction to forcing file:
            //std::println("   seaice_min = {}", seaice_min);
            //std::println("   seaice_max = {}", seaice_max);
            v_seaice.putVar({t,0,0}, {1,ny,nx}, seaice_curr.data());
            // Correct fluxes:
            for (int j=0; j<ny; ++j)
                for (int i=0; i<nx; ++i)
                    albedo[j][i] = (seaice_curr[j][i]*albedo_ice)+((1.0-seaice_curr[j][i])*albedo_ocean);
            v_swrad_frc.getVar({t,0,0}, {1,ny,nx}, swrad_frc.data());
            for (int j=0; j<ny; ++j)
                for (int i=0; i<nx; ++i) {
                    tmp = swrad_frc[j][i]/(1-albedo[j][i]);  // swrad_in in Matlab code
                    swrad_frc[j][i] = tmp*(1-seaice_curr[j][i])*(1.0f-albedo_ocean);
                    swrad_ice_frc[j][i] = tmp*seaice_curr[j][i]*(1.0f-albedo_ice);
                }
            v_swflux_frc.getVar({t,0,0}, {1,ny,nx}, swflux_frc.data());
            for (int j=0; j<ny; ++j)
                for (int i=0; i<nx; ++i)
                    swflux_frc[j][i] += 1.3*swflux_seaice_curr[j][i];
            // Heat flux:
            v_shflux_frc.getVar({t,0,0}, {1,ny,nx}, shflux_frc.data());
            for (int j=0; j<ny; ++j)
                for (int i=0; i<nx; ++i) {
                    tmp = shflux_frc[j][i];  // shflux_atm in Matlab code
                    shflux_ice_frc[j][i] = tmp*seaice_curr[j][i];  // heat flux into sea ice
                    shflux_frc[j][i] = tmp*(1.0-seaice_curr[j][i]);  // heat flux that goes directly into the ocean
                }
            for (int j=0; j<ny; ++j)
                for (int i=0; i<nx; ++i) {
                    shflux_ice_lat = 1.3*shflux_seaice_curr[j][i]; // latent heat flux associated with melting and freezing
                    shflux_ice_melt = shflux_ice_lat>0 ? 0.0f : shflux_ice_lat;
                    //shflux_ice_melt(shflux_ice_lat>0)=0;
                    //shflux_ice_freeze=shflux_ice_lat;
                    //shflux_ice_freeze(shflux_ice_lat<0)=0;
                    shflux_ice_freeze = shflux_ice_lat<0 ? 0.0f : shflux_ice_lat;
                    // Ice melts from below and freezes from above:
                    // i.e. we add the latent heat of freezing to the total ice
                    // flux and also the short wave flux on top of the ice. 
                    // The net is sensible heat into the ice the acts to
                    // cool or warm the ice layer (ice heat capacity). This heat
                    // stored in the ice is released upon melting to the ocean (will
                    // cool the ocean if negative).
                    shflux_ice_frc[j][i] += shflux_ice_freeze + swrad_ice_frc[j][i];
                    shflux_melt = shflux_ice_melt + shflux_snowice_curr[j][i];  // total latent heat from sea ice and snow melt
                    // Scale 90% of sensible flux (stored heat) with total melt, 10% heat conduction:
                    shflux_ice_sen = 0.9 * shflux_melt * nansum_shflux_ice/nansum_shflux_melt + 0.1*shflux_ice_frc[j][i];
                    //shflux_ice_sen=(shflux_melt./...
                    //    repmat(nansum(shflux_melt(:)),size(shflux_ice)).*...
                    //    repmat(nansum(shflux_ice(:).*0.9),size(shflux_ice)))+shflux_ice.*0.1; 
                    // Add latent heat of melt and sensible heat flux under ice:
                    shflux_frc[j][i] += shflux_ice_melt + shflux_ice_sen;
                }
            // Wind stress:
            // Check if sustr is on the u or rho grid:
            v_sustr = nc_frc.getVar("sustr");
            nx_sustr = nx;
            sustr_on_ugrid = false;
            for (auto &dim : v_sustr.getDims()) {
                if (dim.getName() == "xi_u") sustr_on_ugrid = true;
            }
            if (sustr_on_ugrid) nx_sustr = nx - 1;
            array_2d sustr(boost::extents[ny][nx_sustr]);
            v_sustr.getVar({t,0,0}, {1,ny,nx_sustr}, sustr.data());
            if (sustr_on_ugrid) {
                for (int j=0; j<ny; ++j)
                    for (int i=0; i<nx_sustr; ++i) {
                        tmp1 = 0.5f * (seaice_curr[j][i]+seaice_curr[j][i+1]);
                        tmp2 = sustr[j][i] * tmp1 * (1.0f - tmp1);
                        sustr[j][i] = sustr[j][i]*(1.0-tmp1) + tmp2;
                    }
            } else {
                for (int j=0; j<ny; ++j)
                    for (int i=0; i<nx; ++i) {
                        tmp = sustr[j][i] * seaice_curr[j][i] * (1.0f - seaice_curr[j][i]);
                        sustr[j][i] = sustr[j][i]*(1.0-seaice_curr[j][i]) + tmp;
                    }
            }
            v_sustr.putVar({t,0,0}, {1,ny,nx_sustr}, sustr.data());
            // Check if svstr is on the v or rho grid:
            v_svstr = nc_frc.getVar("svstr");
            ny_svstr = ny;
            svstr_on_vgrid = false;
            for (auto &dim : v_svstr.getDims()) {
                if (dim.getName() == "eta_v") svstr_on_vgrid = true;
            }
            if (svstr_on_vgrid) ny_svstr = ny - 1;
            array_2d svstr(boost::extents[ny_svstr][nx]);
            v_svstr.getVar({t,0,0}, {1,ny_svstr,nx}, svstr.data());
            if (svstr_on_vgrid) {
                for (int j=0; j<ny_svstr; ++j)
                    for (int i=0; i<nx; ++i) {
                        tmp1 = 0.5f * (seaice_curr[j][i]+seaice_curr[j+1][i]);
                        tmp2 = svstr[j][i] * tmp1 * (1.0f - tmp1);
                        svstr[j][i] = svstr[j][i]*(1.0-tmp1) + tmp2;
                    }
            } else {
                for (int j=0; j<ny; ++j)
                    for (int i=0; i<nx; ++i) {
                        tmp = svstr[j][i] * seaice_curr[j][i] * (1.0f - seaice_curr[j][i]);
                        svstr[j][i] = svstr[j][i]*(1.0-seaice_curr[j][i]) + tmp;
                    }
            }
            v_svstr.putVar({t,0,0}, {1,ny_svstr,nx}, svstr.data());
            // Finalize shflux and swflux:
            for (int j=0; j<ny; ++j)
                for (int i=0; i<nx; ++i) {
                    swflux_frc[j][i] += swflux_snowice_curr[j][i];
                    shflux_frc[j][i] += shflux_snowice_curr[j][i];
                }
            // Write corrected fluxes to output file:
            v_swflux_frc.putVar({t,0,0}, {1,ny,nx}, swflux_frc.data());
            v_shflux_frc.putVar({t,0,0}, {1,ny,nx}, shflux_frc.data());
            v_swrad_frc.putVar({t,0,0}, {1,ny,nx}, swrad_frc.data());
        }
        // Add some variable attributes to forcing file:
        attval = seaice_file.c_str();
        v_swflux_frc.putAtt("seaice_corr", NcType::nc_CHAR, strlen(attval), attval);
        v_shflux_frc.putAtt("seaice_corr", NcType::nc_CHAR, strlen(attval), attval);
        v_swrad_frc.putAtt("seaice_corr", NcType::nc_CHAR, strlen(attval), attval);
        v_seaice.putAtt("seaice_corr", NcType::nc_CHAR, strlen(attval), attval);
        attval = snowice_file.c_str();
        v_swflux_frc.putAtt("snowice_corr", NcType::nc_CHAR, strlen(attval), attval);
        v_shflux_frc.putAtt("snowice_corr", NcType::nc_CHAR, strlen(attval), attval);
        // Add global attribute indicating that the sea ice correction has been
        // done on this file:
        attval = "completed";
        nc_frc.putAtt("ROMSpy_adjustment_seaice_correction", NcType::nc_CHAR, strlen(attval), attval);
        nc_frc.close();
    }

    std::println("------- END of ROMS_forcing::make_seaice_correction -------\n");
    std::cout.flush();
}