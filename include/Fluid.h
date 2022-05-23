#ifndef FLUID_H
#define FLUID_H

#include <iostream>
#include <cmath>
#include <Eigen/Dense>
#include <Eigen/Sparse>
#include<Eigen/IterativeLinearSolvers>
#include <vector>
#include <Eigen/StdVector>
#include <iterator>
#include "ppties.h"
#include "Mesh.h"

using namespace Eigen;

class Fluid
{
	MatrixXf vsol, vtemp;
	properties fluid_ppts;
	VectorXf vcelerity, num_elems, vpres;

	vector<float> presL2t;

	Mesh fl_mesh, fl_mesh_np1;

public:
	Fluid(properties ppts_);
	float timestep(VectorXf vcor, MatrixXf vsol, VectorXf wx, VectorXf vcelerity, float CFL);
	void initialize(Mesh msh);
	void solve(float Delta_t, float u_dot_t);
	void Lax_Wendroff(float gamm1, float Delta_t, VectorXf wx, Mesh msh_n, Mesh msh_np1);
	MatrixXf flu_residual(int ie, float Delta_t, VectorXf vcore0, MatrixXi conec, MatrixXf vflue, VectorXf wxe);
	MatrixXf flux(MatrixXf vsol, VectorXf vpres, VectorXf wx);
	MatrixXf shock_capture(VectorXf vcor, MatrixXi conec, MatrixXf vmgn, MatrixXf vmgnp1,
						   MatrixXf &vres, MatrixXf vsol, VectorXf xlumpm, VectorXf &num_elems);
	void store_data(vector<float> &histo_pressure, vector<float> &Imp_fl, vector<VectorXf, aligned_allocator<VectorXf> > &histo_velocity,
					vector<VectorXf, aligned_allocator<VectorXf> > &histo_pres_field,
					vector<VectorXf, aligned_allocator<VectorXf> > &histo_rho, vector<VectorXf, aligned_allocator<VectorXf> > &histo_rho_v,
					vector<VectorXf, aligned_allocator<VectorXf> > &histo_rho_e, vector<VectorXf, aligned_allocator<VectorXf> > &histo_mesh, 
					Mesh &msh, float Delta_t, float u_dot_t, int istep);
	MatrixXf flu_mass(int ie, VectorXf vcore);
	VectorXf get_vpres();
	MatrixXf get_vsol();
	VectorXf get_vcelerity();

	friend class FSI;
};

#endif
