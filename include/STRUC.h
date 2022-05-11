#ifndef STRUC_H
#define STRUC_H

#include <iostream>
#include <cmath>
#include <Eigen/Dense>
#include <vector>
#include <Eigen/StdVector>
#include <iterator>
#include "ppties.h"
#if defined(_LINUX) | (_WIN32)
#include <pybind11/pybind11.h>
#include <pybind11/embed.h>
#include <pybind11/eigen.h>
#endif

using namespace Eigen;
#if defined(_LINUX) | (_WIN32)
namespace py = pybind11;
using namespace py::literals;
#endif
class STRUC
{

    float u_t, u_dot_t, u_double_dot_t, Ppiston, delta_u;
    properties struc_ppts;

public:
    float freq0, omega0, T0, dt_export;
    #if defined(_LINUX) | (_WIN32)
    py::object drom;
    #endif
    STRUC(properties ppt);
    void set_ppts(properties ppt);
    void set_BC(float presL2t_ind);
    void solve(float Delta_T);
    void lin_model_solve(float Delta_t);
    void nonlin_model_solve(float Delta_t);
    void rom_model_solve(float Delta_t);
    void store_data(vector<VectorXf, aligned_allocator<VectorXf> > &histo_deformation,
                    vector<float> &histo_accel, vector<float> &Force_ext, vector<float> &Ec,
                    vector<float> &Ep, vector<float> &Em);
    float get_u();
    float get_u_dot_t();
    float get_u_double_dot_t();
    float get_Ppiston();
    void initialize(float presPist);

    friend class FSI;
};

#endif
