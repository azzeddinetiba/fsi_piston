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
#include <Eigen/StdVector>
#include <fstream>
#include <iterator>


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
				 float &e_init, MatrixXf &vsol, VectorXf &vpres, VectorXf &vtemp, VectorXf &vcelerity, float &cel_init, 
				 vector<float> &presL2t)
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
	presL2t.push_back(pres_init);

}

void data_s_init(int &nnt, float &Lsp0, float &pres_init0, float &Lspe, float &p_ext,
				 float &A, float &omega0, float &freq0, float &T0, VectorXf &vpres,
				 float &pres_init, float &vsols0, float &u_t, vector<float> &vprel, float &u_double_dot_t)
{

	float vfg0;

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

void spring(vector<float> vprel, float Ppiston, float pres_init0, float A, float Delta_t,
			float vsols0, float &u_t, float &u_dot_t, float &u_double_dot_t)
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

void structure(int &nnt, float &vsols0, int &istep, vector<float> &presL2t,
			   float &u_t, float &u_dot_t, float &u_double_dot_t, vector<float> &vprel,
			   float &pres_init0, float &A, float &Delta_t, float &Lspe, float &Lsp0,
			   vector<VectorXf, aligned_allocator<VectorXf> > &histo_deformation, vector<float> &Ec, vector<float> &Ep,
			   vector<float> &Em, vector<float> &Force_ext, VectorXf &vpres)
{
	vsols0 = u_t;
	float Ppiston;

	if (istep == 0)
	{
		Ppiston = presL2t[0];
	}
	else
	{
		Ppiston = presL2t[istep - 1];
	}

	spring(vprel, Ppiston, pres_init0, A, Delta_t, vsols0, u_t, u_dot_t, u_double_dot_t);

	Ec.push_back(.5 * vprel[1] * pow(u_dot_t, 2));
	Ep.push_back(.5 * vprel[0] * pow((Lspe - u_t - Lsp0), 2));
	Em.push_back(Ec[istep] + Ep[istep]);

	Force_ext.push_back(-(vpres(nnt - 1) - vpres(0)) * A);

	VectorXf hist_def(2);
	hist_def(0) = u_t;
	hist_def(1) = 0;
	histo_deformation.push_back(hist_def);
}

void move_mesh(float &Delta_t, float u_dot_t, MatrixXf vsol,
			   VectorXf &vcor_n, VectorXf &vcor_np1, VectorXf &wx)
{
	int nnt, ndim;
	nnt = vcor_n.rows();
	ndim = vcor_n.cols();
	wx = (u_dot_t * VectorXf::LinSpaced(Sequential, nnt, 0, 1));

	vcor_n = vcor_np1;

	vcor_np1.col(0) = vcor_n.col(0) + Delta_t * wx;
}

MatrixXf flux(MatrixXf vsol, VectorXf vpres, VectorXf wx)
{
	ArrayXXf foo(vsol.rows(), 3), vsol_;
	ArrayXf vpres_, wx_;
	vsol_ = vsol.array();
	vpres_ = vpres.array();
	wx_ = wx.array();

	foo.col(0) = vsol_.col(1) - vsol_.col(0) * wx_;
	foo.col(1) = vsol_.col(1).pow(2) / vsol_.col(0) + vpres_ - vsol_.col(1) * wx_;
	foo.col(2) = (vsol_.col(2) + vpres_) * vsol_.col(1) / vsol_.col(0) - vsol_.col(2) * wx_;

	return foo.matrix();
}

MatrixXf flu_residual(int &ie, float &Delta_t, MatrixXi &conec, MatrixXf &vsol,
					  VectorXf &vpres, VectorXf &vcore0, MatrixXf &vflue, VectorXf &wxe,
					  float &gamm1)
{
	VectorXi kloce;
	MatrixXf vsole, vsnod(2, vsol.cols()), vflu1_2, vrese(2, 3), interm(2, 2), interm2(2, 3);
	VectorXf vprese, x, vsol1_2, vsolmoy, vpres1_2SV;
	ArrayXf vpres1_2S;
	ArrayXXf vsnod_;
	float xl0;

	kloce = conec(ie, all);
	vsole = vsol(kloce, all);
	vprese = vpres(kloce);

	x = vcore0(1, all) - vcore0(0, all);
	xl0 = sqrt(x.transpose() * x);

	vsol1_2 = .5 * (vsole(0, all) + vsole(1, all)) -
			  .5 * Delta_t / xl0 * (vflue(1, all) - vflue(0, all));

	vsolmoy = (vsole(0, all) + vsole(1, all)) * .5;

	vsnod = vsole;
	vsnod.row(0) += vsol1_2.transpose() - vsolmoy.transpose();
	vsnod.row(1) += vsol1_2.transpose() - vsolmoy.transpose();
	vsnod_ = vsnod.array();

	vpres1_2S = gamm1 * (vsnod_(all, 2) - .5 * vsnod_(all, 1).pow(2) / vsnod_(all, 0));
	vpres1_2SV = vpres1_2S.matrix();

	vflu1_2 = flux(vsnod, vpres1_2SV, wxe);

	interm(0, 0) = -1;
	interm(0, 1) = -1;
	interm(1, 0) = 1;
	interm(1, 1) = 1;

	interm2.row(0) = vflu1_2.row(0);
	interm2.row(1) = vflu1_2.row(1);
	vrese = .5 * Delta_t * (interm * interm2);

	return vrese;
}

MatrixXf flu_mass(int &ie, VectorXf &vcore)
{
	float xl;
	VectorXf x;
	MatrixXf vme, interm(2, 2);

	x = vcore(1, all) - vcore(0, all);

	xl = sqrt(x.transpose() * x);
	interm(0, 0) = 2;
	interm(0, 1) = 1;
	interm(1, 0) = 1;
	interm(1, 1) = 2;
	vme = xl / 6 * interm;

	return vme;
}

MatrixXf shock_capture(VectorXf vcor, MatrixXi &conec, MatrixXf &vmgn, MatrixXf &vmgnp1,
					   MatrixXf &vres, MatrixXf &vsol, VectorXf &xlumpm, VectorXf &num_elems)
{
	
	int nnt, ndln, nelt, nnel, num_el;
	float cd, xl;
	nnt = vsol.rows();
	ndln = vsol.cols();
	nelt = conec.rows();
	nnel = conec.cols();

	VectorXf pmax(nnt), pmin(nnt), rmax(nnt), rmin(nnt), vresj, vsolj, vdul, vduh, 
		uimax, uimin, x, vcore, qmax, qmin, var_cel(nnt), celi(2);
	MatrixXf vduh3(nnt, 3), vdul3(nnt, 3), aec(nelt, 2), interm(2, 2), cel(nelt, 2), cel1(nelt, 2), 
		cel2(nelt, 2), cel3(nelt, 2), du(nnt, ndln);
	VectorXi kloce;
	LLT<MatrixXf> llt_vmgnp1;
	cel1 = MatrixXf::Zero(nnt, 2);
	cel2 = MatrixXf::Zero(nnt, 2);
	cel3 = MatrixXf::Zero(nnt, 2);
	pmax = VectorXf::Zero(nnt);
	pmin = VectorXf::Zero(nnt);
	rmax = VectorXf::Zero(nnt);
	rmin = VectorXf::Zero(nnt);

	interm(0, 0) = 1;
	interm(0, 1) = -1;
	interm(1, 0) = -1;
	interm(1, 1) = 1;

	cd = .4;

	llt_vmgnp1.compute(vmgnp1);


	for (int j = 0; j < ndln; j++)
	{
		vresj = vres.col(j);
		vsolj = vsol.col(j);

		vresj = vresj + (vmgn - vmgnp1) * vsolj;

		vdul = ((vresj + cd * vmgn * vsolj - cd * (xlumpm.array() * vsolj.array()).matrix()).array() / xlumpm.array()).matrix();
		vduh = llt_vmgnp1.solve(vresj);

		
		vduh3.col(j) = vduh;
		vdul3.col(j) = vdul;


		uimax = -1e10 * VectorXf::Ones(nnt);
		uimin = 1e10 * VectorXf::Ones(nnt);


		for (int ie = 0; ie < nelt; ie++)
		{
			kloce = conec(ie, all);
			vcore = vcor(kloce);
			x = vcore(1, all) - vcore(0, all);
			xl = sqrt(x.transpose() * x);

			uimax(kloce) = uimax(kloce).cwiseMax((vsolj(kloce).cwiseMax(vsolj(kloce) + vdul(kloce))).maxCoeff());
			uimin(kloce) = uimin(kloce).cwiseMin((vsolj(kloce).cwiseMin(vsolj(kloce) + vdul(kloce))).minCoeff());

			aec(ie, all) = ((xl / 6 * interm * (cd * vsolj(kloce) + vduh(kloce))).array() / xlumpm(kloce).array()).matrix();
			pmax(kloce) = pmax(kloce) + (aec(ie, all)).cwiseMax(VectorXf::Zero(2).transpose()).transpose();
			pmin(kloce) = pmin(kloce) + (aec(ie, all)).cwiseMin(VectorXf::Zero(2).transpose()).transpose();
		}

		qmax = uimax - (vsolj + vdul);
		qmin = uimin - (vsolj + vdul);


		for (int in = 0; in < nnt; in++)
		{
			if ( (std::abs(pmin(in)) < 1e-9) || (std::abs(pmax(in)) < 1e-9))
			{
				rmax(in) = 0;
				rmin(in) = 0;
			}
			else
			{
				rmax(in) = min(1.f, qmax(in) / pmax(in));
				rmin(in) = min(1.f, qmin(in) / pmin(in));
			}
		}



			
		for (int ie = 0; ie < nelt; ie++)
		{
			for (int in = 0; in < 2; in++)
			{


				if (aec(ie, in) >= 0)
				{
					cel(ie, in) = min(rmax(conec(ie, 0)), rmax(conec(ie, 1)));
				}
				else
				{
					cel(ie, in) = min(rmin(conec(ie, 0)), rmin(conec(ie, 1)));
				}

			}
		}


		if (j == 0) cel1 = cel;

		if (j == 1) cel2 = cel;

		if (j == 2) cel3 = cel;

	}
				cout<<"\n";


	var_cel = VectorXf::Zero(nnt);
	for (int ie = 0; ie < nelt; ie++)
	{
		celi = (cel1(ie, all).cwiseMin(cel2(ie, all))).cwiseMin(cel3(ie, all));
		var_cel(conec(ie, all)) = var_cel(conec(ie, all)) + (celi.array() / num_elems(conec(ie, all)).array()).matrix();
	}

	for (int j = 0; j < ndln; j++)
	{
		du(all, j) = vdul3(all, j) + (var_cel.array() * (vduh3(all, j) - vdul3(all, j)).array()).matrix();
	}


	return du;
}

void Lax_Wendroff(float &gamm1, float &Delta_t, MatrixXf &vsol, VectorXf &vpres,
				  VectorXf &wx, MatrixXi &conec, VectorXf &vcor_n, VectorXf &vcor_np1, VectorXf &num_elems)
{
	int nnt, ndln, nelt, nnel, ie;
	nnt = vsol.rows();
	nelt = conec.rows();
	ndln = vsol.cols();
	nnel = conec.cols();

	MatrixXf vres = MatrixXf::Zero(nnt, ndln), vmg_n = MatrixXf::Zero(nnt, nnt),
			 vmg_np1 = MatrixXf::Zero(nnt, nnt), vflux, vflue, vrese, du(nnt, ndln);
	VectorXi kloce;
	VectorXf vcore_n, vcore_np1, wxe, xlumpm(nnt);

	vflux = flux(vsol, vpres, wx);

	for (ie = 0; ie < nelt; ie++)
	{
		kloce = conec(ie, seq(0, nnel - 1));
		vcore_n = vcor_n(kloce.array(), 0);
		vcore_np1 = vcor_np1(kloce.array(), 0);
		vflue = vflux(kloce.array(), all);
		wxe = wx(kloce);

		vrese = flu_residual(ie, Delta_t, conec, vsol, vpres,
							 vcore_n, vflue, wxe, gamm1);


		vres(kloce, all) = vres(kloce, all) + vrese;

		vmg_n(kloce, kloce) = vmg_n(kloce, kloce) + flu_mass(ie, vcore_n);

		vmg_np1(kloce, kloce) = vmg_np1(kloce, kloce) + flu_mass(ie, vcore_np1);

		if (ie == 0)
		{
			vres.row(0) = vres.row(0) + Delta_t * vflue.row(0);
		}

		if (ie == nelt - 1)
		{
			vres.row(nnt - 1) = vres.row(nnt - 1) - Delta_t * vflue.row(1);
		}
	}


	xlumpm = vmg_n.colwise().sum().transpose();

	du = shock_capture(vcor_n, conec, vmg_n, vmg_np1, vres, vsol, xlumpm, num_elems);

	du(0, 1) = 0;
	du(nnt - 1, 1) = 0;

	vsol = vsol + du;
}

void fluid(int &nnt, float &Delta_t, float &u_dot_t, MatrixXf &vsol, VectorXf &vcor_n,
		   VectorXf &vcor_np1, VectorXf &wx, float &gamm1, VectorXf &vpres, MatrixXi &conec,
		   VectorXf &num_elems, VectorXf &vtemp, float &C_v, float &R, int &istep,
		   vector<float> &presL2t, float &A, vector<float> &Imp_fl,
		   vector<VectorXf, aligned_allocator<VectorXf> > &histo_velocity, vector<float> &histo_pressure,
		   vector<VectorXf, aligned_allocator<VectorXf> > &histo_deformation, float &L_0)
{

	VectorXf hist_veloc(2), wx_vec(1);
	move_mesh(Delta_t, u_dot_t, vsol, vcor_n, vcor_np1, wx);

	vsol(nnt - 1, 1) = vsol(nnt - 1, 0) * u_dot_t;

	Lax_Wendroff(gamm1, Delta_t, vsol, vpres, wx, conec, vcor_n, vcor_np1, num_elems);

	vtemp = temperature(C_v, vsol);
	vpres = pressure(R, vtemp, vsol);
	presL2t.push_back(vpres(nnt - 1));

	if (istep == 0)
		Imp_fl.push_back(presL2t[istep] * A * u_dot_t * Delta_t);
	else
		Imp_fl.push_back(Imp_fl[istep - 1] + presL2t[istep] * A * u_dot_t * Delta_t);


	hist_veloc(0) = u_dot_t;
	hist_veloc(1) = wx(nnt - 1);
	histo_velocity.push_back(hist_veloc);
	histo_pressure.push_back(presL2t[istep]);
	histo_deformation[istep](1) = vcor_np1(nnt - 1, 0) - L_0;
}

int main()
{

	// Physical properties
	float gam, gamm1, R, C_v, pres_init0, temp_init0, rho_init0, pres_init, rho_init, temp_init, p_ext, u_init, e_init;
	float L_0, L_t, U_0, A;
	int nelt, nnt, nnel, ndln, ndle, ndlt;
	float cel_init, Lsp0, Lspe, omega0, freq0, T0, vsols0, u_t = .2;
	nelt = 70;
	nnt = nelt + 1;
	nnel = 2;
	ndln = 3;
	ndle = ndln * nnel;
	ndlt = nnt * ndln;
	float u_dot_t = 0., vfg0, u_double_dot_t;

	// Connectivity
	MatrixXi conec(nelt, 2);
	conec.col(0) = VectorXi::LinSpaced(Sequential, nelt, 0, nnt - 2);
	conec.col(1) = VectorXi::LinSpaced(Sequential, nelt, 1, nnt - 1);

	// Initialisation
	MatrixXf vflux(nnt, ndln);
	vector<float> Delta_t_storage, t;
	vector<float> Ec, Ep, Em, Force_ext, presL2t;
	VectorXf x, vcor, vcor0, vcor_n, vcor_np1;
	VectorXf vpres, vtemp, vcelerity;
	MatrixXf vsol(nnt, 3);

	VectorXf wx(nnt);
	vector<float> vprel, histo_pressure, Imp_fl;
	VectorXf temp, pres;

	vector<VectorXf, aligned_allocator<VectorXf> > histo_velocity,
		histo_deformation;

	data_f_init(nelt, nnt, nnel, ndln, ndle, ndlt, L_0,
				A, U_0, L_t, vcor, vcor0, vcor_n, vcor_np1,
				gam, gamm1, R, C_v, pres_init0, temp_init0, rho_init0,
				pres_init, rho_init, temp_init, p_ext, u_init,
				e_init, vsol, vpres, vtemp, vcelerity, cel_init, presL2t);
	data_s_init(nnt, Lsp0, pres_init0, Lspe, p_ext, A, omega0, freq0,
				T0, vpres, pres_init, vsols0, u_t, vprel, u_double_dot_t);

	float Tmax, CFL, dxmin, Delta_t;
	Tmax = .5 * T0;
	CFL = .8;

	dxmin = (vcor_np1(seq(1, nnt - 1)).array() - vcor_np1(seq(0, nnt - 2)).array()).minCoeff();
	Delta_t = CFL * dxmin / (vcelerity.array() + (vsol.array().col(1) / vsol.array().col(0)).abs()).maxCoeff();

	VectorXf num_elems(nnt);
	for (int ind = 0; ind < nnt; ind++)
	{
		if (ind == 0 || ind == nnt - 1)
			num_elems(ind) = 1;
		else
			num_elems(ind) = 2;
	}

	float Total_time = 0;
	int istep = -1;

	while (Total_time < (Tmax - Delta_t))
	{
		istep += 1;
		Delta_t = timestep(vcor_np1, vsol, wx, vcelerity, CFL);
		Delta_t_storage.push_back(Delta_t);
		Total_time += Delta_t;
		t.push_back(Total_time);

		cout<<" \n Total Time \n";
		cout<<Total_time;

		structure(nnt, vsols0, istep, presL2t,
				  u_t, u_dot_t, u_double_dot_t, vprel, pres_init0, A, Delta_t,
				  Lspe, Lsp0, histo_deformation, Ec, Ep, Em, Force_ext, vpres);
		

		fluid(nnt, Delta_t, u_dot_t, vsol, vcor_n, vcor_np1, wx, gamm1, vpres, conec, num_elems, 
			vtemp, C_v, R, istep, presL2t, A, Imp_fl, histo_velocity, histo_pressure, 
			histo_deformation, L_0);
	}

	cout << "\n Succeeded";

	std::ofstream file("results_sol.txt");
	if (file.is_open())
	{
		file << "Here is the matrix solution:\n" << vsol << '\n';
	}

	std::ofstream file2("results_pres.txt");
	if (file2.is_open())
	{
		    std::ostream_iterator<float> output_iterator(file2, "\n");
			std::copy(histo_pressure.begin(), histo_pressure.end(), output_iterator);
	}
	return 0;
}
