/*
 * main.cpp
 *
 *  Created on: Dec 10, 2021
 *      Author: tiba
 */

//

#define _USE_MATH_DEFINES

#include <iostream>
#include <cmath>
#include <Eigen/Dense>
#include <vector>
#include <Eigen/StdVector>
#include <fstream>
#include <iterator>
#include "STRUC.h"
#include "Mesh.h"
#include "Fluid.h"
#include "FSI.h"
#include "config.h"
#include <pybind11/pybind11.h>
#include <pybind11/embed.h>
#include <pybind11/eigen.h>

using namespace Eigen;
using namespace std;

properties load_ppts()
{
	properties ppts;
	float a, b, c, interm;

	ppts.Coeff = coeff; // Fraction of natural structural period, giving the total period of simulation

	ppts.L_0 = 1; // Initial Gas Chamber Length (not including the initial displacement)
	ppts.A = 1;	  // Section

	ppts.U_0 = .2; // Initial displacement
	ppts.L_t = ppts.L_0 + ppts.U_0;

	ppts.gam = 1.4; // the specific heat ratio of the gas
	ppts.gamm1 = ppts.gam - 1.;
	ppts.R = 287;					// the individual gas constant
	ppts.C_v = ppts.R / ppts.gamm1; // the specific heat capacity of the gas

	ppts.pres_init0 = 1E5;														// initial pressure for chamber length = L0
	ppts.temp_init0 = 300;														// initial temperature
	ppts.rho_init0 = ppts.pres_init0 / ppts.gamm1 / ppts.C_v / ppts.temp_init0; // initial volumic mass

	ppts.pres_init = ppts.pres_init0 * pow((ppts.L_0 / ppts.L_t), ppts.gam);

	ppts.rho_init = ppts.rho_init0 * pow((ppts.pres_init / ppts.pres_init0), (1. / ppts.gam));
	ppts.temp_init = ppts.pres_init / ppts.rho_init / ppts.gamm1 / ppts.C_v;
	ppts.p_ext = 0 * ppts.pres_init0; // pressure on the right of the piston

	// we set the initial fluid velocity and the initial total fluid energy
	ppts.u_init = 0.;
	ppts.e_init = ppts.pres_init / ppts.gamm1 / ppts.rho_init + 0.5 * pow(ppts.u_init, 2.);

	ppts.vprel.push_back(1e7);	// Spring rigidity
	ppts.vprel.push_back(mass); // Spring mass
	ppts.sdim = 1;
	ppts.spring_model = "linear";
	ppts.nln_order = 3;
	ppts.rom_in_struc = false;
	ppts.cont_rom = true;

	ppts.Lsp0 = 1.2;								   // Unstretched spring length
	ppts.young = ppts.vprel[0] * ppts.Lsp0 / ppts.A;   // Young modulus equivalent to the spring rigidity
	ppts.rho_s = ppts.vprel[1] / (ppts.A * ppts.Lsp0); // Beam density equivalent to the spring mass
	if (ppts.spring_model == "nonlinear")
	{
		ppts.umax = 0.2; // Maximum spring displacements for linear spring model ('C' Model)
		ppts.mu = mu_coeff * ppts.vprel[0] / ppts.umax;
		if (ppts.nln_order == 2)
		{
			ppts.u0 = (-ppts.vprel[0] + sqrt(pow(ppts.vprel[0], 2) + 4 * ppts.mu * ppts.A * ppts.pres_init0)) / (-2 * ppts.mu);
		}
		else
		{
			a = ppts.vprel[0];
			b = ppts.mu;
			c = ppts.A * ppts.pres_init0;
			interm = pow((((std::sqrt((27 * b * pow(c, 2) + 4 * pow(a, 3)) / b)) / (b * 2 * pow(3, (3. / 2.)))) - c / (2 * b)), (1. / 3.));
			ppts.u0 = interm - a / (3 * b * interm);
		}
		ppts.Lspe = ppts.Lsp0 + ppts.u0;
	}
	else
	{
		ppts.Lspe = ppts.Lsp0 - (ppts.pres_init0 - ppts.p_ext) * ppts.A / ppts.vprel[0]; // initial spring length
	}

	ppts.dt = 7.09e-6;

	// Newmark params
	ppts.newm = true;
	ppts.newm_gamma = .5;
	ppts.newm_beta = .25 * std::pow(ppts.newm_gamma + .5, 2);

	// Generalized alpha scheme
	ppts.ch_alph = false;
	ppts.ch_rho = 0.;
	ppts.ch_alpha_m = (2 * ppts.ch_rho - 1.) / (ppts.ch_rho + 1.);
	ppts.ch_alpha_f = ppts.ch_rho / (ppts.ch_rho + 1);
	ppts.ch_beta = .25 * std::pow(1. - ppts.ch_alpha_m + ppts.ch_alpha_f, 2);
	ppts.ch_gamma = .5 - ppts.ch_alpha_m + ppts.ch_alpha_f;

	return ppts;
}

int main()
{

	// Geometrical and physical properties
	properties ppts;
	ppts = load_ppts();
	float dt = 1e-5, Tmax = 0.18, Total_time = 0., Delta_t;
	vector<float> t;
	vector<VectorXf, aligned_allocator<VectorXf>> histo_un;
	vector<VectorXf, aligned_allocator<VectorXf>> histo_udt;
	vector<VectorXf, aligned_allocator<VectorXf>> histo_uddt;
	ofstream file1("../../test_results/results_un.txt");
	ofstream file2("../../test_results/results_udt.txt");
	ofstream file3("../../test_results/results_uddt.txt");
	VectorXf vcor, vcor_np1, vcelerity, wx;
	MatrixXf vsol;
	float u_dot_t, ppiston;

	// Create the meshes
	int nnt = nmesh, i = 0;
	Mesh mesh_ns;
	mesh_ns.load(nnt, ppts.Lsp0);

	// Create the structure FEM model
	STRUC structure_model(ppts);
	structure_model.initialize(0., mesh_ns);
	structure_model.store_test_data(histo_un, histo_udt, histo_uddt);
	t.push_back(Total_time);

	if (file1.is_open())
	{
		file1 << histo_un[0] << '\n';
	}

	if (file2.is_open())
	{
		file2 << histo_udt[0] << '\n';
	}

	if (file3.is_open())
	{
		file3 << histo_uddt[0] << '\n';
	}

	Delta_t = dt;

	// Time loop/increments
	while (Total_time < (Tmax - Delta_t))
	{

		cout << "Total Time :\n";
		cout << Total_time;
		cout << "\n";

		// Solve the structural problem

		// Apply the piston pressure value from the fluid problem to the solid problem
		structure_model.set_BC(0.);

		// Get the structure solution
		structure_model.solve(Delta_t);
		structure_model.store_test_data(histo_un, histo_udt, histo_uddt);

		if (file1.is_open())
		{
			file1 << histo_un[i+1] << '\n';
		}

		if (file2.is_open())
		{
			file2 << histo_udt[i+1] << '\n';
		}

		if (file3.is_open())
		{
			file3 << histo_uddt[i+1] << '\n';
		}
		// Storing the time data
		Total_time += Delta_t;
		t.push_back(Total_time);

		i += 1;
	}

	ofstream file("../../test_results/results_t.txt");
	if (file.is_open())
	{
		std::ostream_iterator<float> output_iterator(file, "\n");
		std::copy(t.begin(), t.end(), output_iterator);
	}

	return 0;
}