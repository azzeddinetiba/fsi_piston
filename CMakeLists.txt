cmake_minimum_required(VERSION 3.16.0)
project(fsi_piston)

#Parameters
set (coeff 2.5 CACHE STRING "")
set (mass 1000. CACHE STRING "")
set (nmesh 250 CACHE STRING "")
set (mu_coeff 6 CACHE STRING "")
configure_file (include/config.h.in ../include/config.h @ONLY)

#Include libraries
include_directories(include)
find_package (Eigen3 3.3 REQUIRED NO_MODULE)

#Add source files
file(GLOB SOURCES "src/*.cpp")

add_executable(fsi_piston ${SOURCES})
target_link_libraries (fsi_piston Eigen3::Eigen)
