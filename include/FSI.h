#ifndef FSI_H
#define FSI_H

#include <iostream>
#include <cmath>
#include <Eigen/Dense>
#include <vector>
#include <Eigen/StdVector>
#include <iterator>
#include <fstream>
#include "ppties.h"
#include "Mesh.h"
#include "STRUC.h"
#include "Fluid.h"

using namespace Eigen;

class FSI
{
	public:
		vector<float> histo_pressure, Force_ext, Ec, Ep, Em;
		vector<float> Imp_fl;
		vector<VectorXf, aligned_allocator<VectorXf> > histo_velocity;
		vector<VectorXf, aligned_allocator<VectorXf> > histo_deformation;
		vector<VectorXf, aligned_allocator<VectorXf> > histo_pres_field;
		vector<VectorXf, aligned_allocator<VectorXf> > histo_rho;
		vector<VectorXf, aligned_allocator<VectorXf> > histo_rho_v;
		vector<VectorXf, aligned_allocator<VectorXf> > histo_rho_e;

		float Tmax, CFL, dxmin, Delta_t, Total_time;
		vector<float> t, Delta_t_storage;
		int istep;

		FSI(float T0);
		void solve(STRUC &struc, Fluid &fluid);
		void export_results();


};

#endif
