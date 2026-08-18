#include <nanobind/nanobind.h>
#include <nanobind/stl/string.h>
#include <nanobind/stl/vector.h>
#include "roms_forcing_utils.h"

namespace nb = nanobind;
using namespace nb::literals;

// nanobind code for Python binding:
NB_MODULE(roms_forcing_utils, m) {
    nb::class_<ROMS_forcing>(m, "ROMS_frc")
        .def(nb::init<const std::string,const std::string,const std::string,bool>())
        .def("rel_humidity", &ROMS_forcing::rel_humidity)
        .def("wind_stress", &ROMS_forcing::wind_stress)
        .def("make_river_freshwater", &ROMS_forcing::make_river_freshwater)
        .def("make_seaice_correction", &ROMS_forcing::make_seaice_correction);
    //m.def("rel_humidity", &rel_humidity, "outfile"_a, "year"_a, "t2m_file"_a, "d2m_file"_a,
    //    "Computes relative humidity according to the recommendation of the ERA5 documentation.");
    m.doc() = "Class for computing certain parts of a ROMS forcing";
}