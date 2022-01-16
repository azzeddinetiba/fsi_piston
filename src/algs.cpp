#define _USE_MATH_DEFINES

#include "algs.h"

newton::newton(int iter_max, float tol)
{
    max_iter = iter_max;
    tolerance = tol;
}

void newton::initialize(float init)
{
    iter = 0;
    dx = 1e4;
    initial_delta_x = init;
}

bool newton::criterion(float delta_x)
{
    bool criter;
    if (iter == 0)
        criter = true;
    else
        criter = iter < max_iter && (abs(dx) / abs(delta_x)) > tolerance;

    if (!criter && iter != 0)
    {
        std::cout << "\nSystem converged at iteration " << iter << ", residual = " << residual << std::endl;
    }
    return criter;
}

float newton::advance()
{
    dx = -residual / jacobian;
    iter += 1;
    return -residual / jacobian;
}

void newton::set_jacobian(float jac)
{
    jacobian = jac;
}

void newton::set_residual(float res)
{
    residual = res;
}

float newton::get_init()
{
    return initial_delta_x;
}
