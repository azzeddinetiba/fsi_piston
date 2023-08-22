cmake_minimum_required(VERSION 3.16.0)
project(wave_equ_tst)

#Parameters
set (coeff 2.5 CACHE STRING "")
set (mass 1000. CACHE STRING "")
set (nmesh 250 CACHE STRING "")
set (mu_coeff 6 CACHE STRING "")
configure_file (../include/config.h.in ../../include/config.h @ONLY)

#Include libraries
include_directories(../include)
find_package (Eigen3 3.3 REQUIRED NO_MODULE)
if (CMAKE_SYSTEM_NAME STREQUAL "Linux")
    find_package(pybind11 REQUIRED)
endif ()
if (CMAKE_SYSTEM_NAME STREQUAL "Windows")
    add_subdirectory(pybind11)
    link_directories("C:/Users/azzed/anaconda3/libs")
endif ()

#Add source files
add_executable(wave_equ_tst ../src/algs.cpp ../src/Mesh.cpp ../src/STRUC.cpp test_wave.cpp)

if (CMAKE_SYSTEM_NAME STREQUAL "Darwin")
    target_link_libraries(wave_equ_tst PRIVATE Eigen3::Eigen)
else()
    target_link_libraries(wave_equ_tst PRIVATE Eigen3::Eigen pybind11::pybind11)
endif ()
