#include <netcdf>
#include <vector>
#include <cmath>
#include <boost/multi_array.hpp>

#pragma once

using namespace netCDF;

class ROMS_forcing {
public:
    ROMS_forcing(std::string,std::string,std::string,bool);
    void rel_humidity(std::string outfile, int year, std::string t2m_file, std::string d2m_file);
    void wind_stress(std::string outfile, int year, std::string era5_ewss_file, std::string era5_nsss_file,
        std::string ref_file, std::string time_dim_ref, int tstep_day);
    std::string make_river_freshwater(std::string river_input_file, std::string outdir, int time_res,
        std::string spreading);
    void save_river_freshwater(std::string out_file, std::string tfreq, boost::multi_array<float,3> data,
        std::string river_input_file, int days_per_year);
    void make_seaice_correction(std::vector<std::string> frc_files, int year, std::string seaice_file,
        std::string snowice_file, int time_res_h, bool verbose);
private:
    std::string ROMS_setup;
    std::string ROMS_gridfile;
    std::string output_dir;
    bool verbose = true;
};