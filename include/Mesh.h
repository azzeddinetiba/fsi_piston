#ifndef MESH_H
#define MESH_H

#include <iostream>
#include <cmath>
#include <Eigen/Dense>
#include <vector>
#include <Eigen/StdVector>
#include <iterator>
#include "ppties.h"

using namespace Eigen;


class Mesh
{

	VectorXf vcor, wx;
	MatrixXi conec;

public:
	int nnt, nelt;
	void load(int nnt_, float L_t);
	void move_mesh(float Delta_t, float u_dot_t, MatrixXf vsol);
	VectorXf get_vcor();
	VectorXf get_wx();
	MatrixXi get_conec();

};

#endif
