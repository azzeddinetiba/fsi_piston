#define _USE_MATH_DEFINES

#include "FSI.h"

using namespace Eigen;
using namespace std;

FSI::FSI(float T0)
{
	Total_time = 0;
	istep = -1;
	CFL = .8;
}

void FSI::export_results()
{

	ofstream file("../Results/results_pres.txt");
	if (file.is_open())
	{
		std::ostream_iterator<float> output_iterator(file, "\n");
		std::copy(histo_pressure.begin(), histo_pressure.end(), output_iterator);
	}

	ofstream file6("../Results/results_Delta_T.txt");
	if (file6.is_open())
	{
		std::ostream_iterator<float> output_iterator(file6, "\n");
		std::copy(Delta_t_storage.begin(), Delta_t_storage.end(), output_iterator);
	}

	ofstream file1("../Results/results_pres_field.txt");
	ofstream file2("../Results/results_rho.txt");
	ofstream file3("../Results/results_rho_v.txt");
	ofstream file4("../Results/results_rho_e.txt");
	ofstream file5("../Results/results_v.txt");
	ofstream file7("../Results/results_mesh.txt");
	ofstream file8("../Results/results_m_accel.txt");
	ofstream file9("../Results/results_Ec.txt");
	ofstream file10("../Results/results_Ep.txt");
	ofstream file11("../Results/results_Em.txt");
	ofstream file12("../Results/results_Imp_Fl.txt");

	for (int i = 0; i < istep + 1; i++)
	{
		if (file1.is_open())
		{
			file1 << histo_pres_field[i] << '\n';
		}

		if (file2.is_open())
		{
			file2 << histo_rho[i] << '\n';
		}

		if (file3.is_open())
		{
			file3 << histo_rho_v[i] << '\n';
		}

		if (file4.is_open())
		{
			file4 << histo_rho_e[i] << '\n';
		}

		if (file5.is_open())
		{
			file5 << histo_velocity[i] << '\n';
		}

		if (file7.is_open())
		{
			file7 << histo_mesh[i] << '\n';
		}

		if (file8.is_open())
		{
			file8 << histo_accel[i] << '\n';
		}
		if (file9.is_open())
		{
			file9 << Ec[i] << '\n';
		}
		if (file10.is_open())
		{
			file10 << Ep[i] << '\n';
		}
		if (file11.is_open())
		{
			file11 << Em[i] << '\n';
		}
		if (file12.is_open())
		{
			file12 << Imp_fl[i] << '\n';
		}
	}
}

void FSI::solve(STRUC &struc, Fluid &fluid, float d_t)
{

	VectorXf vcor, vcor_np1, vcelerity, wx;
	MatrixXf vsol;
	int nnt;
	float u_dot_t, ppiston;

	vcor_np1 = fluid.fl_mesh_np1.get_vcor();
	nnt = fluid.fl_mesh_np1.nnt;
	vsol = fluid.get_vsol();
	Tmax = fluid.fluid_ppts.Coeff * struc.T0;

	vcelerity = fluid.get_vcelerity();
	dxmin = (vcor_np1(seq(1, nnt - 1)).array() - vcor_np1(seq(0, nnt - 2)).array()).minCoeff();
	if (d_t < 1e-9)
		Delta_t = CFL * dxmin / (vcelerity.array() + (vsol.array().col(1) / vsol.array().col(0)).abs()).maxCoeff();
	else
		Delta_t = d_t;

	// Time loop/increments
	while (Total_time < (Tmax - Delta_t))
	{

		cout << "Total Time :\n";
		cout << Total_time;
		cout << "\n";

		istep += 1;

		// Declare the fluid solution (pressure at the piston)
		vcor_np1 = fluid.fl_mesh_np1.get_vcor();
		wx = fluid.fl_mesh_np1.get_wx();
		vsol = fluid.get_vsol();

		// Compute the next time step value respecting CFL
		if (d_t < 1e-9)
			Delta_t = fluid.timestep(vcor_np1, vsol, wx, vcelerity, CFL);
		else
			Delta_t = d_t;
		Delta_t_storage.push_back(Delta_t);
		Total_time += Delta_t;

		t.push_back(Total_time);

		// Solve the structural problem

		// Apply the piston pressure value from the fluid problem to the solid problem
		if (istep == 0)
		{
			ppiston = fluid.presL2t[0];
		}
		else
		{
			ppiston = fluid.presL2t[istep - 1];
		}
		struc.set_BC(ppiston);

		// Get the structure solution
		struc.solve(Delta_t);
		u_dot_t = struc.get_u_dot_t();
		struc.store_data(histo_deformation, histo_accel, Force_ext, Ec, Ep, Em);

		// Solve the fluid problem
		fluid.solve(Delta_t, u_dot_t);
		fluid.store_data(histo_pressure, Imp_fl, histo_velocity,
						 histo_deformation, histo_pres_field, histo_rho, histo_rho_v,
						 histo_rho_e, histo_mesh, fluid.fl_mesh_np1, Delta_t, u_dot_t, istep);
	}
}
