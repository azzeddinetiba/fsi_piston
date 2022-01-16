#ifndef algs_H
#define algs_H

#include <iostream>
#include <cmath>

class newton
{
    int iter, max_iter;
    float tolerance, jacobian, residual, dx, initial_delta_x;

public:
    newton(int iter_max = 500, float tol = 1e-5);
    bool criterion(float delta_x);
    float advance();
    void initialize(float init = 0.);
    void set_residual(float res);
    void set_jacobian(float jac);
    float get_init();
};

#endif
