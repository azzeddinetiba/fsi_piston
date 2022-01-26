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
	if (struc_ppts.spring_model == "nonlinear")
	{
		Ep.push_back(.5 * struc_ppts.vprel[0] * pow((u0 - u_t), 2));

		if ((u0 - u_t) > 0)
			Ep[Ep.size() - 1] += mu * (1 / 3) * pow((u0 - u_t), 3);
		else
			Ep[Ep.size() - 1] -= mu * (1 / 3) * pow((u0 - u_t), 3);
	}
	else
	{
		Ep.push_back(.5 * struc_ppts.vprel[0] * pow((struc_ppts.Lspe - u_t - struc_ppts.Lsp0), 2));
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
	float vkg, vmg, vfg, vkt, vres, inter, A;
	float res, jacobian, du;
	int i = 0;
	newton nonlin_sys;

	A = struc_ppts.A;

	nonlin_sys.initialize();
	delta_u = nonlin_sys.get_init();
	std::cout << "\n Newton iterations ... ";
	std::cout << "\n";
	while (nonlin_sys.criterion(delta_u))
	{

		res = pow(Delta_t, 2) * (A * Ppiston + struc_ppts.vprel[0] * u0 + struc_ppts.vprel[1] * u_double_dot_t -
								 struc_ppts.vprel[0] * u_t - struc_ppts.vprel[0] * delta_u - mu * (u_t - u0 + delta_u) * std::abs(u_t + delta_u - u0)) -
			  4 * struc_ppts.vprel[1] * delta_u + 4 * struc_ppts.vprel[1] * u_dot_t * Delta_t;

		if ((u_t - u0 + delta_u) > 0.)
		{
			jacobian = -4 * struc_ppts.vprel[1] + pow(Delta_t, 2) * (-struc_ppts.vprel[0] + 2 * mu * (u_t - u0) + 2 * mu * delta_u);
		}
		else
		{
			jacobian = -4 * struc_ppts.vprel[1] + pow(Delta_t, 2) * (-struc_ppts.vprel[0] - 2 * mu * (u_t - u0) - 2 * mu * delta_u);
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
	// Initialisation of the acceleration
	u_t = struc_ppts.U_0;
	float vfg0;
	vfg0 = (presPist - 0 * struc_ppts.pres_init0) * struc_ppts.A;
	u_dot_t = 0;
	if (struc_ppts.spring_model == "linear")
	{
		u_double_dot_t = (vfg0 + struc_ppts.vprel[0] * (struc_ppts.Lspe - u_t - struc_ppts.Lsp0)) / struc_ppts.vprel[1];
	}
	else
	{
		mu = struc_ppts.vprel[0] / struc_ppts.umax;
		u0 = (-struc_ppts.vprel[0] + sqrt(pow(struc_ppts.vprel[0], 2) + 4 * mu * struc_ppts.A * struc_ppts.pres_init0)) / (-2 * mu);

		u_double_dot_t = (struc_ppts.A * presPist + struc_ppts.vprel[0] * (u0 - u_t) -
						  mu * (u_t - u0) * std::abs(u_t - u0)) /
						 struc_ppts.vprel[1];
	}
}
