#define _USE_MATH_DEFINES

#include "STRUC.h"

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

void STRUC::store_data(vector<VectorXf, aligned_allocator<VectorXf>> &histo_deformation,
					   vector<float> &histo_accel, vector<float> &Force_ext, vector<float> &Ec,
					   vector<float> &Ep, vector<float> &Em)
{
	Ec.push_back(.5 * struc_ppts.vprel[1] * pow(u_dot_t, 2));
	Ep.push_back(.5 * struc_ppts.vprel[0] * pow((struc_ppts.Lspe - u_t - struc_ppts.Lsp0), 2));
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

void STRUC::solve(float Delta_t)
{
	float pres_init0, vkg, vmg, vfg, vkt, vres, vdu, inter, A;
	pres_init0 = struc_ppts.pres_init0;
	A = struc_ppts.A;
	vkg = struc_ppts.vprel[0];
	vmg = struc_ppts.vprel[1];
	vfg = (Ppiston - pres_init0) * A;

	vkt = 4 * vmg + pow(Delta_t, 2) * vkg;

	vres = pow(Delta_t, 2) * vfg + 4 * Delta_t * vmg * u_dot_t +
		   pow(Delta_t, 2) * vmg * u_double_dot_t - pow(Delta_t, 2) * vkg * u_t;

	vdu = vres / vkt;

	u_t += vdu;
	inter = 4 / pow(Delta_t, 2) * vdu - 4 / Delta_t * u_dot_t -
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
	u_double_dot_t = (vfg0 + struc_ppts.vprel[0] * (struc_ppts.Lspe - u_t - struc_ppts.Lsp0)) / struc_ppts.vprel[1];
}
