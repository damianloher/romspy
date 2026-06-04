# Compilation: g++-15 must be available (e.g. on begonia)
# - activate conda env romspy:
#    source activate romspy
# - execute "cmake -B ..." and "cmake --build"

cmake_minimum_required(VERSION 3.15...3.28)

set(CMAKE_CXX_COMPILER "/usr/bin/g++-15")
set(CMAKE_CXX_STANDARD 23)
set(CMAKE_CXX_STANDARD_REQUIRED True)

project(roms_forcing_utils CXX)

# Search for the Development.Module component required by nanobind. The name of
# this module changed across CMake versions, hence the additional conditional check.
if (CMAKE_VERSION VERSION_LESS 3.18)
  set(DEV_MODULE Development)
else()
  set(DEV_MODULE Development.Module)
endif()
find_package(Python 3.8 COMPONENTS Interpreter ${DEV_MODULE} REQUIRED)

# Perform an optimized release build by default unless another build type is specified:
if (NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
  set(CMAKE_BUILD_TYPE Release CACHE STRING "Choose the type of build." FORCE)
  set_property(CACHE CMAKE_BUILD_TYPE PROPERTY STRINGS "Debug" "Release" "MinSizeRel" "RelWithDebInfo")
endif()

# Detect the installed nanobind package and import it into CMake
execute_process(
  COMMAND "${Python_EXECUTABLE}" -m nanobind --cmake_dir
  OUTPUT_STRIP_TRAILING_WHITESPACE OUTPUT_VARIABLE nanobind_ROOT)
find_package(nanobind CONFIG REQUIRED)

# Define source code files of extension:
nanobind_add_module(roms_forcing_utils 
   src/roms_forcing_utils.h
   src/roms_forcing_utils.cpp 
   src/river_forcing.cpp
   src/seaice.cpp
   src/extrapolate_land.cpp
   src/nanobind.cpp)

# The following paths are from module netcdf/4.7.4:
set(NETCDF_ROOT /usr/local/netcdf-4.7.4-4.5.3-4.3.1-gnu-7.5.0)
#set(NETCDF_ROOT /usr/local/netcdf-4.7.3-4.5.2-4.3.1-gnu-7.5.0)
set(NETCDF_INC_PATH ${NETCDF_ROOT}/include)
set(NETCDF_LIB_PATH ${NETCDF_ROOT}/lib64)
target_include_directories(roms_forcing_utils PRIVATE ${NETCDF_INC_PATH})
target_include_directories(roms_forcing_utils PRIVATE nanoflann/include)
target_link_directories(roms_forcing_utils PRIVATE ${NETCDF_LIB_PATH})
target_link_libraries(roms_forcing_utils PRIVATE netcdf_c++4 netcdf)