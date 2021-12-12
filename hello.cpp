/*
 * hello.cpp
 *
 *  Created on: Dec 10, 2021
 *      Author: tiba
 */

// Your First C++ Program

#define _USE_MATH_DEFINES

#include <iostream>
#include <cmath>
#include <Eigen/Dense>
#include <vector>

using namespace Eigen;
using namespace std;

VectorXf temperature(float C_v, MatrixXf vsol)
{

	VectorXf temp(vsol.rows());

	ArrayXf col0, col1, col2, res;
	col0 = vsol(all, 0).array();
	col1 = vsol(all, 1).array();
	col2 = vsol(all, 2).array();

	res = (col2 / col0 - (.5 * col1.square()) / col0.square()) / C_v;

	temp = res.matrix();

	return temp;
}

VectorXf pressure(float R, VectorXf vtemp, MatrixXf vsol)
{

	VectorXf pres(vsol.rows());

	ArrayXf col0, res;
	col0 = vsol(all, 0).array();

	res = col0 * R * vtemp.array();

	pres = res.matrix();

	return pres;
}

void data_f_init(int &nelt, int &nnt, int &nnel, int &ndln, int &ndle, int &ndlt, float &L_0,
				 float &A, float &U_0, float &L_t, VectorXf &vcor, VectorXf &vcor0, VectorXf &vcor_n, VectorXf &vcor_np1,
				 float &gam, float &gamm1, float &R, float &C_v, float &pres_init0, float &temp_init0, float &rho_init0,
				 float &pres_init, float &rho_init, float &temp_init, float &p_ext, float &u_init,
				 float &e_init, MatrixXf &vsol, VectorXf &vpres, VectorXf &vtemp, VectorXf &vcelerity, float &cel_init)
{

	L_0 = 1;
	A = 1;

	U_0 = .2;
	L_t = L_0 + U_0;

	// Vector of coordinates

	vcor = VectorXf::LinSpaced(Sequential, nnt, 0.0, L_t);
	vcor0 = vcor;
	vcor_n = vcor;
	vcor_np1 = vcor;

	// Connectivity
	MatrixXd connec(nelt, 2);
	connec.col(0) = VectorXd::LinSpaced(Sequential, nelt, 0, nnt - 2);
	connec.col(1) = VectorXd::LinSpaced(Sequential, nelt, 1, nnt - 1);

	gam = 1.4; // the specific heat ratio of the gas
	gamm1 = gam - 1.;
	R = 287;		 // the individual gas constant
	C_v = R / gamm1; // the specific heat capacity of the gas

	pres_init0 = 1E5;								   // initial pressure for chamber length = L0
	temp_init0 = 300;								   // initial temperature
	rho_init0 = pres_init0 / gamm1 / C_v / temp_init0; // initial volumic mass

	pres_init = pres_init0 * pow((L_0 / L_t), gam);

	rho_init = rho_init0 * pow((pres_init / pres_init0), (1. / gam));
	temp_init = pres_init / rho_init / gamm1 / C_v;
	p_ext = 0 * pres_init0; // pressure on the right of the piston

	// we set the initial fluid velocity and the initial total fluid energy
	u_init = 0.;
	e_init = pres_init / gamm1 / rho_init + 0.5 * pow(u_init, 2.);

	// initialisation of the solution matrix

	ArrayXXf vsol_i(nnt, 3);
	vsol_i.col(0) = rho_init * VectorXf::Ones(nnt).array();
	vsol_i.col(1) = rho_init * u_init * VectorXf::Ones(nnt).array();
	vsol_i.col(2) = rho_init * e_init * VectorXf::Ones(nnt).array();

	vsol = vsol_i.matrix();

	ArrayXf vcelerity_i;

	vtemp = temperature(C_v, vsol);
	vpres = pressure(R, vtemp, vsol);
	vcelerity_i = (gam * gamm1 * C_v * vtemp.array()).sqrt();
	vcelerity = vcelerity_i.matrix();
	cel_init = vcelerity(0);
}

void data_s_init(int &nnt, float &Lsp0, float &pres_init0, float &Lspe, float &p_ext,
				 float &A, float &omega0, float &freq0, float &T0, VectorXf &vpres,
				 float &pres_init, float &vsols0, float &u_t, vector<float> &vprel)
{

	float vfg0, u_double_dot_t;

	vprel.push_back(1e7);
	vprel.push_back(100);

	Lsp0 = 1.2;
	Lspe = Lsp0 - (pres_init0 - p_ext) * A / vprel[0];

	omega0 = sqrt(vprel[0] / vprel[1]);
	freq0 = omega0 / (2 * M_PI);
	T0 = 1. / freq0;

	// Initialisation of solution for the structure
	vsols0 = u_t;

	// Initialisation of the acceleration
	vfg0 = (vpres(nnt - 1) - 0 * pres_init) * A;
	u_double_dot_t = (vfg0 + vprel[0] * (Lspe - u_t - Lsp0)) / vprel[1];
}

float timestep(VectorXf &vcor, MatrixXf &vsol, VectorXf &wx, VectorXf &vcelerity, float &CFL)
{

	int nnt;
	float dmin, Delta_t, dxmin;

	nnt = vcor.size();
	dxmin = (vcor(seq(1, nnt - 1)) - vcor(seq(0, nnt - 2))).minCoeff();

	Delta_t = CFL * dxmin / ((vsol.array().col(1) / vsol.array().col(0)).abs() + vcelerity.array() + wx.array().abs()).maxCoeff();

	return Delta_t;
}

void spring(vector<float> &vprel, float &Ppiston, float &pres_init0, float &A, float &Delta_t,
			float &vsols0, float &u_t, float &u_dot_t, float &u_double_dot_t)
{
	float vkg, vmg, vfg, vkt, vres, vdu, inter;
	vkg = vprel[0];
	vmg = vprel[1];
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

void structure(int &nnt, float &vsols0, int& istep, MatrixXf &presL2t,
			   float &u_t, float &u_dot_t, float &u_double_dot_t, vector<float> &vprel,
			   float &pres_init0, float &A, float &Delta_t, float &Lspe, float &Lsp0,
			   vector<float> &histo_deformation, vector<float> &Ec, vector<float> &Ep,
			   vector<float> &Em, vector<float> &Force_ext, VectorXf &vpres)
{
	vsols0 = u_t;
	float Ppiston;

	if (istep == 0)
	{
		Ppiston = vpres(nnt - 1);
	}
	else
	{
		Ppiston = presL2t(istep - 1, 2);
	}

	spring(vprel, Ppiston, pres_init0, A, Delta_t, vsols0, u_t, u_dot_t, u_double_dot_t);

	histo_deformation.push_back(u_t);
	Ec.push_back(.5 * vprel[1] * pow(u_dot_t, 2));
	Ep.push_back(.5 * vprel[0] * pow((Lspe - u_t - Lsp0), 2));
	Em.push_back(Ec[istep] + Ep[istep]);

	Force_ext.push_back(-(vpres(nnt - 1) - vpres(0)) * A);
}

void move_mesh(float &Delta_t, float &u_dot_t, MatrixXf &vsol,
			   VectorXf &vcor_n, VectorXf &vcor_np1, VectorXf &wx)
{
	int nnt, ndim;
	nnt = vcor_n.rows();
	ndim = vcor_n.cols();
	wx = (u_dot_t * VectorXf::LinSpaced(Sequential, nnt, 0, 1).array()).matrix();

	vcor_n = vcor_np1;

	vcor_np1.col(0) = vcor_n.col(0) + Delta_t * wx;
}

void Lax_Wendroff(float& gamm1, float& Delta_t, vsol, vpres, wx, conec, vcor_n, vcor_np1, number)
{
	int nnt, ndln, nelt, nnel;
	nnt = vsol.rows(); nelt = conec.rows();
	ndln = vsol.cols(); nnel = conec.cols();

}

void fluid(int& nnt, float &Delta_t, float &u_dot_t, MatrixXf &vsol, VectorXf &vcor_n, 
		VectorXf &vcor_np1, VectorXf &wx, float gamm1, VectorXf vpres, MatrixXf conec)
{

	move_mesh(Delta_t, u_dot_t, vsol, vcor_n, vcor_np1, wx);
	vsol(nnt - 1, 1) = vsol(nnt - 1, 0) * u_dot_t;

	Lax_Wendroff(gamm1, Delta_t, vsol, vpres, wx, conec, vcor_n, vcor_np1, number);


}

int main()
{

	// Physical properties
	float gam, gamm1, R, C_v, pres_init0, temp_init0, rho_init0, pres_init, rho_init, temp_init, p_ext, u_init, e_init;
	float L_0, L_t, U_0, A;
	int nelt, nnt, nnel, ndln, ndle, ndlt;
	float cel_init, Lsp0, Lspe, omega0, freq0, T0, vsols0, u_t = .2;
	nelt = 20;
	nnt = nelt + 1;
	nnel = 2;
	ndln = 3;
	ndle = ndln * nnel;
	ndlt = nnt * ndln;
	float u_dot_t = 0., vfg0, u_double_dot_t;

	// Initialisation
	MatrixXf vflux(nnt, ndln), presL2t;
	vector<float> Delta_t_storage, t;
	vector<float> histo_deformation, Ec, Ep, Em, Force_ext;
	VectorXf x, vcor, vcor0, vcor_n, vcor_np1;
	VectorXf vpres, vtemp, vcelerity;
	MatrixXf vsol(nnt, 3);

	VectorXf wx(nnt);
	vector<float> vprel;
	VectorXf temp, pres;
	MatrixXf m(2, 3);
	m(0, 0) = 3;
	m(1, 0) = 2.5;
	m(0, 1) = -1;

	temp = temperature(1., m);
	pres = pressure(1, temp, m);

	VectorXd x1, x2;
	x1 = VectorXd::LinSpaced(Sequential, 71 - 1, 1, 70);

	data_f_init(nelt, nnt, nnel, ndln, ndle, ndlt, L_0,
				A, U_0, L_t, vcor, vcor0, vcor_n, vcor_np1,
				gam, gamm1, R, C_v, pres_init0, temp_init0, rho_init0,
				pres_init, rho_init, temp_init, p_ext, u_init,
				e_init, vsol, vpres, vtemp, vcelerity, cel_init);
	data_s_init(nnt, Lsp0, pres_init0, Lspe, p_ext, A, omega0, freq0,
				T0, vpres, pres_init, vsols0, u_t, vprel);

	float Tmax, CFL, dxmin, Delta_t;
	Tmax = .5 * T0;
	CFL = .8;

	dxmin = (vcor_np1(seq(1, nnt - 1)).array() - vcor_np1(seq(0, nnt - 2)).array()).minCoeff();
	Delta_t = CFL * dxmin / (vcelerity.array() + (vsol.array().col(1) / vsol.array().col(0)).abs()).maxCoeff();

	float Total_time = 0;
	int istep = -1;

	while (Total_time < (Tmax - Delta_t))
	{
		istep += 1;
		Delta_t = timestep(vcor_np1, vsol, wx, vcelerity, CFL);
		Delta_t_storage.push_back(Delta_t);
		Total_time += Delta_t;
		t.push_back(Total_time);

		structure(nnt, vsols0, istep, presL2t,
				  u_t, u_dot_t, u_double_dot_t, vprel, pres_init0, A, Delta_t,
				  Lspe, Lsp0, histo_deformation, Ec, Ep, Em, Force_ext, vpres);

		fluid(nnt, Delta_t, u_dot_t, vsol, vcor_n, vcor_np1, wx);
	}

	cout<<"\n Succeeded";

	return 0;
}
