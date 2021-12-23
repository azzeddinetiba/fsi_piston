#ifndef STRUC_H
#define STRUC_H

#include <iostream>
#include <cmath>
#include <Eigen/Dense>
#include <vector>
#include <Eigen/StdVector>
#include <iterator>
#include "ppties.h"

using namespace Eigen;

class STRUC
{

    float u_t, u_dot_t, u_double_dot_t, Ppiston;
    properties struc_ppts;

public:
    float freq0, omega0, T0;
    STRUC(properties ppt);
    void set_ppts(properties ppt);
    void set_BC(float presL2t_ind);
    void solve(float Delta_T);
    void store_data(vector<VectorXf, aligned_allocator<VectorXf> > &histo_deformation,
                    vector<float> &Force_ext, vector<float> &Ec, vector<float> &Ep, vector<float> &Em);
    float get_u();
    float get_u_dot_t();
    float get_u_double_dot_t();
    float get_Ppiston();
    void initialize(float presPist);

    friend class FSI;
};

#endif
