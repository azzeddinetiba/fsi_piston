#define _USE_MATH_DEFINES

#include "Fluid.h"

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

VectorXf Fluid::get_vcelerity()
{
	return vcelerity;
}

MatrixXf Fluid::get_vsol()
{
	return vsol;
}

VectorXf Fluid::get_vpres()
{
	return vpres;
}

Fluid::Fluid(properties ppts_)
{
	fluid_ppts = ppts_;
}

void Fluid::initialize(Mesh msh)
{
	fl_mesh = msh;
	fl_mesh_np1 = msh;
	ArrayXXf vsol_i(msh.nnt, 3);
	vsol_i.col(0) = fluid_ppts.rho_init * VectorXf::Ones(msh.nnt).array();
	vsol_i.col(1) = fluid_ppts.rho_init * fluid_ppts.u_init * VectorXf::Ones(msh.nnt).array();
	vsol_i.col(2) = fluid_ppts.rho_init * fluid_ppts.e_init * VectorXf::Ones(msh.nnt).array();

	vsol = vsol_i.matrix();

	presL2t.push_back(fluid_ppts.pres_init);

	vtemp = temperature(fluid_ppts.C_v, vsol);
	vpres = pressure(fluid_ppts.R, vtemp, vsol);

	ArrayXf vcelerity_i;
	vcelerity_i = (fluid_ppts.gam * fluid_ppts.gamm1 * fluid_ppts.C_v * vtemp.array()).sqrt();
	vcelerity = vcelerity_i.matrix();

	num_elems = VectorXf::Zero(msh.nnt);
	for (int ind = 0; ind < msh.nnt; ind++)
	{
		if (ind == 0 || ind == msh.nnt - 1)
			num_elems(ind) = 1;
		else
			num_elems(ind) = 2;
	}
}

float Fluid::timestep(VectorXf vcor, MatrixXf vsol, VectorXf wx, VectorXf vcelerity, float CFL)
{

	int nnt;
	float dmin, Delta_t, dxmin;

	nnt = vcor.size();
	dxmin = (vcor(seq(1, nnt - 1)) - vcor(seq(0, nnt - 2))).minCoeff();

	Delta_t = CFL * dxmin / ((vsol.array().col(1) / vsol.array().col(0)).abs() + vcelerity.array() + wx.array().abs()).maxCoeff();

	return Delta_t;
}

MatrixXf Fluid::flux(MatrixXf vsol, VectorXf vpres, VectorXf wx)
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

MatrixXf Fluid::flu_residual(int ie, float Delta_t, VectorXf vcore0, MatrixXi conec, MatrixXf vflue, VectorXf wxe)
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

	vpres1_2S = fluid_ppts.gamm1 * (vsnod_(all, 2) - .5 * vsnod_(all, 1).pow(2) / vsnod_(all, 0));
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

MatrixXf Fluid::shock_capture(VectorXf vcor, MatrixXi conec, MatrixXf vmgn, MatrixXf vmgnp1,
							  MatrixXf &vres, MatrixXf vsol, VectorXf xlumpm, VectorXf &num_elems)
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
	SparseMatrix<float> s_vmgn(nnt, nnt);
	SparseMatrix<float> s_vmgnp1(nnt, nnt);
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

	// llt_vmgnp1.compute(vmgnp1);
	s_vmgn = vmgn.sparseView();
	s_vmgnp1 = vmgnp1.sparseView();

	ConjugateGradient<SparseMatrix<float> > cg_vmgnp1;
	cg_vmgnp1.compute(s_vmgnp1);


	for (int j = 0; j < ndln; j++)
	{
		vresj = vres.col(j);
		vsolj = vsol.col(j);

		vresj = vresj + (s_vmgn - s_vmgnp1) * vsolj;

		vdul = ((vresj + cd * vmgn * vsolj - cd * (xlumpm.array() * vsolj.array()).matrix()).array() / xlumpm.array()).matrix();
		vduh = cg_vmgnp1.solve(vresj);

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
			if ((std::abs(pmin(in)) < 1e-9) || (std::abs(pmax(in)) < 1e-9))
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

		if (j == 0)
			cel1 = cel;

		if (j == 1)
			cel2 = cel;

		if (j == 2)
			cel3 = cel;
	}

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

MatrixXf Fluid::flu_mass(int ie, VectorXf vcore)
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

void Fluid::Lax_Wendroff(float gamm1, float Delta_t, VectorXf wx, Mesh msh_n, Mesh msh_np1)
{

	MatrixXf vres = MatrixXf::Zero(msh_n.nnt, 3), vmg_n = MatrixXf::Zero(msh_n.nnt, msh_n.nnt),
			 vmg_np1 = MatrixXf::Zero(msh_n.nnt, msh_n.nnt), vflux, vflue, vrese, du(msh_n.nnt, 3);
	VectorXi kloce;
	MatrixXi conec;
	VectorXf vcore_n, vcore_np1, wxe, xlumpm(msh_n.nnt), vcor_n, vcor_np1;
	vcor_n = msh_n.get_vcor();
	vcor_np1 = msh_np1.get_vcor();
	conec = msh_n.get_conec();

	vflux = flux(vsol, vpres, wx);

	for (int ie = 0; ie < msh_n.nelt; ie++)
	{
		kloce = conec(ie, seq(0, 1));
		vcore_n = vcor_n(kloce.array(), 0);
		vcore_np1 = vcor_np1(kloce.array(), 0);
		vflue = vflux(kloce.array(), all);
		wxe = wx(kloce);

		vrese = flu_residual(ie, Delta_t, vcore_n, conec, vflue, wxe);

		vres(kloce, all) = vres(kloce, all) + vrese;

		vmg_n(kloce, kloce) = vmg_n(kloce, kloce) + flu_mass(ie, vcore_n);

		vmg_np1(kloce, kloce) = vmg_np1(kloce, kloce) + flu_mass(ie, vcore_np1);

		if (ie == 0)
		{
			vres.row(0) = vres.row(0) + Delta_t * vflue.row(0);
		}

		if (ie == msh_n.nelt - 1)
		{
			vres.row(msh_n.nnt - 1) = vres.row(msh_n.nnt - 1) - Delta_t * vflue.row(1);
		}
	}

	xlumpm = vmg_n.colwise().sum().transpose();

	du = shock_capture(vcor_n, conec, vmg_n, vmg_np1, vres, vsol, xlumpm, num_elems);

	du(0, 1) = 0;
	du(msh_n.nnt - 1, 1) = 0;

	vsol = vsol + du;
}

void Fluid::solve(float Delta_t, float u_dot_t)
{

	VectorXf hist_veloc(2);

	fl_mesh = fl_mesh_np1;

	fl_mesh_np1.move_mesh(Delta_t, u_dot_t, vsol);

	vsol(fl_mesh.nnt - 1, 1) = vsol(fl_mesh.nnt - 1, 0) * u_dot_t;

	Lax_Wendroff(fluid_ppts.gamm1, Delta_t, fl_mesh_np1.get_wx(), fl_mesh, fl_mesh_np1);

	vtemp = temperature(fluid_ppts.C_v, vsol);
	vpres = pressure(fluid_ppts.R, vtemp, vsol);
	presL2t.push_back(vpres(fl_mesh.nnt - 1));
}

void Fluid::store_data(vector<float> &histo_pressure, vector<float> &Imp_fl, vector<VectorXf, aligned_allocator<VectorXf>> &histo_velocity,
					   vector<VectorXf, aligned_allocator<VectorXf>> &histo_deformation, vector<VectorXf, aligned_allocator<VectorXf>> &histo_pres_field,
					   vector<VectorXf, aligned_allocator<VectorXf>> &histo_rho, vector<VectorXf, aligned_allocator<VectorXf>> &histo_rho_v,
					   vector<VectorXf, aligned_allocator<VectorXf>> &histo_rho_e, Mesh &msh, float Delta_t, float u_dot_t, int istep)
{

	int nnt = msh.nnt;
	VectorXf hist_veloc(2);
	VectorXf wx = msh.get_wx();

	if (istep == 0)
		Imp_fl.push_back(presL2t[0] * fluid_ppts.A * u_dot_t * Delta_t);
	else
		Imp_fl.push_back(Imp_fl[istep - 1] + presL2t[istep] * fluid_ppts.A * u_dot_t * Delta_t);

	hist_veloc(0) = u_dot_t;
	hist_veloc(1) = wx(nnt - 1);

	histo_velocity.push_back(hist_veloc);
	histo_pressure.push_back(presL2t[presL2t.size() - 1]);
	histo_deformation[histo_deformation.size() - 1](1) = msh.get_vcor()(nnt - 1, 0) - fluid_ppts.L_0;
	histo_pres_field.push_back(get_vpres());

	histo_rho.push_back(get_vsol().col(0));
	histo_rho_v.push_back(get_vsol().col(1));
	histo_rho_e.push_back(get_vsol().col(2));
}
