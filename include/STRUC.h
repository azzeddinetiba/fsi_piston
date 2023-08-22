#ifndef STRUC_H
#define STRUC_H

#include <iostream>
#include <cmath>
#include <Eigen/Dense>
#include <vector>
#include <Eigen/StdVector>
#include <Eigen/Sparse>
#include <iterator>
#include "ppties.h"
#include "Mesh.h"
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

    float u_t, u_dot_t, u_double_dot_t, Ppiston, delta_u, newm_beta, newm_gamma, ch_alpha_m, ch_alpha_f,
        ch_beta, ch_gamma;
    properties struc_ppts;
    VectorXf rhs, u_n, u_dt_n, u_ddt_n, delta_u_n, u_int_n;
    Eigen::SimplicialCholesky<SparseMatrix<float>> chol_G;
    bool is_there_cholG;
    SparseMatrix<float> rigid, mass;
    Mesh msh;

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
    void lin_1D_model_solve(float Delta_t, float beta, float gamma);
    void nonlin_model_solve(float Delta_t);
    void rom_model_solve(float Delta_t);
    void store_data(vector<VectorXf, aligned_allocator<VectorXf>> &histo_deformation,
                    vector<float> &histo_accel, vector<float> &Force_ext, vector<float> &Ec,
                    vector<float> &Ep, vector<float> &Em, vector<VectorXf, aligned_allocator<VectorXf>> &histo_udt,
                    vector<VectorXf, aligned_allocator<VectorXf>> &histo_uddt);
    void store_test_data(vector<VectorXf, aligned_allocator<VectorXf>> &histo_deformation,
                         vector<VectorXf, aligned_allocator<VectorXf>> &histo_dt,
                         vector<VectorXf, aligned_allocator<VectorXf>> &histo_ddt);
    float get_u();
    float get_u_dot_t();
    float get_u_double_dot_t();
    float get_Ppiston();
    void initialize(float presPist, Mesh mesh);
    void set_BC_essential(int id);
    MatrixXf rigid_e(VectorXf x);
    MatrixXf mass_e(VectorXf x);
    void rhs_term(float p);
    void assemble();

    friend class FSI;
};

#endif
