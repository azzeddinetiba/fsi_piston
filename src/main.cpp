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
	ppts.spring_model = "linear";
	ppts.nln_order = 3;

	ppts.Lsp0 = 1.2; // Unstretched spring length
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

	return ppts;
}

int main()
{

	// Geometrical and physical properties
	properties ppts;
	ppts = load_ppts();

	// Create the mesh
	int nnt = nmesh;
	Mesh mesh_n;
	mesh_n.load(nnt, ppts.L_t);

	// Create the fluid FEM model
	Fluid fluid_model(ppts);
	fluid_model.initialize(mesh_n);

	// Create the structure FEM model
	STRUC structure_model(ppts);
	structure_model.initialize(fluid_model.get_vpres()(nnt - 1));

	// Create the fluid-strucure interaction coupling
	FSI fsi_piston(structure_model.T0);

	// Solve the problem
	fsi_piston.solve(structure_model, fluid_model);

	// Export the results into .txt files
	fsi_piston.export_results();

	return 0;
}
