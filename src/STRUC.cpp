#define _USE_MATH_DEFINES

#include "STRUC.h"
#include "algs.h"

using namespace Eigen;
using namespace std;

STRUC::STRUC(properties ppt)
{
	set_ppts(ppt);
	omega0 = sqrt(struc_ppts.vprel[0] / struc_ppts.vprel[1]);
	freq0 = omega0 / (2 * M_PI);
	T0 = 1. / freq0;
}

float STRUC::get_u()
{
	return u_t;
}

float STRUC::get_u_dot_t()
{
	return u_dot_t;
}

float STRUC::get_u_double_dot_t()
{
	return u_double_dot_t;
}

float STRUC::get_Ppiston()
{
	return Ppiston;
}

void STRUC::store_data(vector<VectorXf, aligned_allocator<VectorXf> > &histo_deformation,
					   vector<float> &histo_accel, vector<float> &Force_ext, vector<float> &Ec,
					   vector<float> &Ep, vector<float> &Em)
{
	Ec.push_back(.5 * struc_ppts.vprel[1] * pow(u_dot_t, 2));
	Ep.push_back(.5 * struc_ppts.vprel[0] * pow((struc_ppts.Lspe - u_t - struc_ppts.Lsp0), 2));
	if (struc_ppts.spring_model == "nonlinear")
	{
		if (struc_ppts.nln_order == 2)
		{
			if ((struc_ppts.u0 - u_t) > 0)
				Ep[Ep.size() - 1] += struc_ppts.mu * (1 / 3) * pow((struc_ppts.u0 - u_t), 3);
			else
				Ep[Ep.size() - 1] -= struc_ppts.mu * (1 / 3) * pow((struc_ppts.u0 - u_t), 3);
		}
		else
		{
			Ep[Ep.size() - 1] += struc_ppts.mu * pow(struc_ppts.u0 - u_t, 4) / 4;
		}
	}

	Em.push_back(Ec[Ec.size() - 1] + Ep[Ep.size() - 1]);
	histo_accel.push_back(u_double_dot_t);

	VectorXf hist_def(2);
	hist_def(0) = u_t;
	hist_def(1) = 0;
	histo_deformation.push_back(hist_def);
}

void STRUC::set_BC(float presL2t_ind)
{
	Ppiston = presL2t_ind;
}

void STRUC::set_ppts(properties ppt)
{
	struc_ppts = ppt;
}

void STRUC::lin_model_solve(float Delta_t)
{
	float pres_init0, vkg, vmg, vfg, vkt, vres, inter, A;
	pres_init0 = struc_ppts.pres_init0;
	A = struc_ppts.A;
	vkg = struc_ppts.vprel[0];
	vmg = struc_ppts.vprel[1];
	vfg = (Ppiston - pres_init0) * A;

	vkt = 4 * vmg + pow(Delta_t, 2) * vkg;

	vres = pow(Delta_t, 2) * vfg + 4 * Delta_t * vmg * u_dot_t +
		   pow(Delta_t, 2) * vmg * u_double_dot_t - pow(Delta_t, 2) * vkg * u_t;

	delta_u = vres / vkt;
}

void STRUC::nonlin_model_solve(float Delta_t)
{
	float vkg, vmg, vfg, vkt, vres, inter, A, mu, kp, u0, mass;
	float res, jacobian, du;
	int i = 0;
	newton nonlin_sys;

	A = struc_ppts.A;
	mu = struc_ppts.mu;
	kp = struc_ppts.vprel[0];
	u0 = struc_ppts.u0;
	mass = struc_ppts.vprel[1];

	nonlin_sys.initialize();
	delta_u = nonlin_sys.get_init();
	std::cout << "\n Newton iterations ... ";
	std::cout << "\n";
	while (nonlin_sys.criterion(delta_u))
	{

		if (struc_ppts.nln_order == 2)
		{
			res = pow(Delta_t, 2) * (A * Ppiston + kp * u0 + mass * u_double_dot_t -
									 kp * u_t - kp * delta_u - mu * (u_t - u0 + delta_u) * std::abs(u_t + delta_u - u0)) -
				  4 * mass * delta_u + 4 * mass * u_dot_t * Delta_t;

			if ((u_t - struc_ppts.u0 + delta_u) > 0.)
			{
				jacobian = -4 * mass + pow(Delta_t, 2) * (-kp + 2 * mu * (u_t - u0) + 2 * mu * delta_u);
			}
			else
			{
				jacobian = -4 * mass + pow(Delta_t, 2) * (-kp - 2 * mu * (u_t - u0) - 2 * mu * delta_u);
			}
		}
		else
		{
			res = mu * pow(delta_u, 3) + (3 * mu * u_t - 3 * mu * u0) * pow(delta_u, 2) +
				  (4 * mass / (pow(Delta_t, 2)) + kp + 3 * mu * pow(u_t, 2) - 6 * mu * u0 * u_t + 3 * mu * pow(u_t, 2)) * delta_u -
				  mass * u_double_dot_t - 4 * mass / (Delta_t)*u_dot_t + (kp + 3 * mu * pow(u0, 2)) * u_t + mu * pow(u_t, 3) -
				  3 * mu * u0 * pow(u_t, 2) - mu * pow(u0, 3) - kp * u0 - A * Ppiston;

			jacobian = 3 * mu * pow(delta_u, 2) + 2 * (3 * mu * u_t - 3 * mu * u0) * delta_u + 4 * mass / (pow(Delta_t, 2)) +
					   kp + 3 * mu * pow(u_t, 2) - 6 * mu * u0 * u_t + 3 * mu * pow(u0, 2);
		}

		nonlin_sys.set_residual(res);
		nonlin_sys.set_jacobian(jacobian);
		delta_u += nonlin_sys.advance();
	}
}

void STRUC::solve(float Delta_t)
{
	float inter;
	if (struc_ppts.spring_model == "linear")
	{
		lin_model_solve(Delta_t);
	}
	else
	{
		nonlin_model_solve(Delta_t);
	}

	u_t += delta_u;
	inter = 4 / pow(Delta_t, 2) * delta_u - 4 / Delta_t * u_dot_t -
			u_double_dot_t;
	u_dot_t = u_dot_t + Delta_t / 2 * (u_double_dot_t + inter);
	u_double_dot_t = inter;
}

void STRUC::initialize(float presPist)
{
	// Initialisation of the displacement
	u_t = struc_ppts.U_0;
	float vfg0;
	vfg0 = (presPist - 0 * struc_ppts.pres_init0) * struc_ppts.A;

	// Initialisation of the velocity
	u_dot_t = 0;

	// Initialisation of the acceleration
	if (struc_ppts.spring_model == "linear")
	{
		u_double_dot_t = (vfg0 + struc_ppts.vprel[0] * (struc_ppts.Lspe - u_t - struc_ppts.Lsp0)) / struc_ppts.vprel[1];
	}
	else
	{
		if (struc_ppts.nln_order == 2)
		{
			u_double_dot_t = (struc_ppts.A * presPist + struc_ppts.vprel[0] * (struc_ppts.u0 - u_t) -
							  struc_ppts.mu * (u_t - struc_ppts.u0) * std::abs(u_t - struc_ppts.u0)) /
							 struc_ppts.vprel[1];
		}
		else
		{
			u_double_dot_t = (struc_ppts.A * presPist + struc_ppts.vprel[0] * (struc_ppts.u0 - u_t) + struc_ppts.mu * pow(struc_ppts.u0 - u_t, 3)) /
							 struc_ppts.vprel[1];
		}
	}
}
