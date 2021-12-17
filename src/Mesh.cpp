#define _USE_MATH_DEFINES

#include "Mesh.h"

void Mesh::load(int nnt_, float L_t)
{
	nnt = nnt_;
	nelt = nnt - 1;
	MatrixXi conec_(nelt, 2);
	conec_.col(0) = VectorXi::LinSpaced(Sequential, nelt, 0, nnt - 2);
	conec_.col(1) = VectorXi::LinSpaced(Sequential, nelt, 1, nnt - 1);
	vcor = VectorXf::LinSpaced(Sequential, nnt, 0.0, L_t);
	conec = conec_;
	wx = VectorXf::Zero(nnt);
}

void Mesh::move_mesh(float Delta_t, float u_dot_t, MatrixXf vsol)

{
	wx = (u_dot_t * VectorXf::LinSpaced(Sequential, nnt, 0, 1));
	vcor.col(0) = vcor.col(0) + Delta_t * wx;
}

VectorXf Mesh::get_vcor()
{
	return vcor;
}

VectorXf Mesh::get_wx()

{
	return wx;
}

MatrixXi Mesh::get_conec()
{
	return conec;
}
