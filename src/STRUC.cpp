#define _USE_MATH_DEFINES

#include "STRUC.h"
#include "algs.h"

using namespace Eigen;

STRUC::STRUC(properties ppt)
{
	set_ppts(ppt);
	if (struc_ppts.sdim == 0)
		omega0 = std::sqrt(struc_ppts.vprel[0] / struc_ppts.vprel[1]);
	else
		omega0 = M_PI * std::sqrt((struc_ppts.young) / (struc_ppts.rho_s)) / (2 * struc_ppts.Lsp0);
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
					   vector<float> &Ep, vector<float> &Em, vector<VectorXf, aligned_allocator<VectorXf>> &histo_dt,
					   vector<VectorXf, aligned_allocator<VectorXf>> &histo_ddt)
{
	Ec.push_back(.5 * struc_ppts.vprel[1] * std::pow(u_dot_t, 2));
	Ep.push_back(.5 * struc_ppts.vprel[0] * std::pow((struc_ppts.Lspe - u_t - struc_ppts.Lsp0), 2));
	if (struc_ppts.spring_model == "nonlinear")
	{
		if (struc_ppts.nln_order == 2)
		{
			if ((struc_ppts.u0 - u_t) > 0)
				Ep[Ep.size() - 1] += struc_ppts.mu * (1 / 3) * std::pow((struc_ppts.u0 - u_t), 3);
			else
				Ep[Ep.size() - 1] -= struc_ppts.mu * (1 / 3) * std::pow((struc_ppts.u0 - u_t), 3);
		}
		else
		{
			Ep[Ep.size() - 1] += struc_ppts.mu * std::pow(struc_ppts.u0 - u_t, 4) / 4;
		}
	}

	Em.push_back(Ec[Ec.size() - 1] + Ep[Ep.size() - 1]);
	histo_accel.push_back(u_double_dot_t);

	histo_deformation.push_back(u_n);
	histo_dt.push_back(u_dt_n);
	histo_ddt.push_back(u_ddt_n);
}

void STRUC::set_BC(float presL2t_ind)
{
	Ppiston = presL2t_ind;
	if (struc_ppts.sdim == 1)
	{
		rhs_term(-Ppiston);
	}
}

void STRUC::set_BC_essential(int id)
{
	for (int i = 0; i < msh.nnt; i++)
	{
		mass.coeffRef(id, i) = 0.;
		mass.coeffRef(i, id) = 0.;
		rigid.coeffRef(id, i) = 0.;
		rigid.coeffRef(i, id) = 0.;
	}
	mass = mass.pruned();
	rigid = rigid.pruned();
	mass.coeffRef(id, id) = 1.;
	rigid.coeffRef(id, id) = 1.;
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

	vkt = 4 * vmg + std::pow(Delta_t, 2) * vkg;

	vres = std::pow(Delta_t, 2) * vfg + 4 * Delta_t * vmg * u_dot_t +
		   std::pow(Delta_t, 2) * vmg * u_double_dot_t - std::pow(Delta_t, 2) * vkg * u_t;

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
			res = std::pow(Delta_t, 2) * (A * Ppiston + kp * u0 + mass * u_double_dot_t -
										  kp * u_t - kp * delta_u - mu * (u_t - u0 + delta_u) * std::abs(u_t + delta_u - u0)) -
				  4 * mass * delta_u + 4 * mass * u_dot_t * Delta_t;

			if ((u_t - struc_ppts.u0 + delta_u) > 0.)
			{
				jacobian = -4 * mass + std::pow(Delta_t, 2) * (-kp + 2 * mu * (u_t - u0) + 2 * mu * delta_u);
			}
			else
			{
				jacobian = -4 * mass + std::pow(Delta_t, 2) * (-kp - 2 * mu * (u_t - u0) - 2 * mu * delta_u);
			}
		}
		else
		{
			res = mu * std::pow(delta_u, 3) + (3 * mu * u_t - 3 * mu * u0) * std::pow(delta_u, 2) +
				  (4 * mass / (std::pow(Delta_t, 2)) + kp + 3 * mu * std::pow(u_t, 2) - 6 * mu * u0 * u_t + 3 * mu * std::pow(u_t, 2)) * delta_u -
				  mass * u_double_dot_t - 4 * mass / (Delta_t)*u_dot_t + (kp + 3 * mu * std::pow(u0, 2)) * u_t + mu * std::pow(u_t, 3) -
				  3 * mu * u0 * std::pow(u_t, 2) - mu * std::pow(u0, 3) - kp * u0 - A * Ppiston;

			jacobian = 3 * mu * std::pow(delta_u, 2) + 2 * (3 * mu * u_t - 3 * mu * u0) * delta_u + 4 * mass / (std::pow(Delta_t, 2)) +
					   kp + 3 * mu * std::pow(u_t, 2) - 6 * mu * u0 * u_t + 3 * mu * std::pow(u0, 2);
		}

		nonlin_sys.set_residual(res);
		nonlin_sys.set_jacobian(jacobian);
		delta_u += nonlin_sys.advance();
	}
}

void STRUC::rom_model_solve(float Delta_t)
{
	VectorXf pred(3), x_input(5);
#if defined(_LINUX) | (_WIN32)
	py::object pred_ob;
#endif

	x_input(0, 0) = u_t;
	x_input(1, 0) = std::pow(u_t, 2);
	x_input(2, 0) = std::pow(u_t, 3);
	x_input(3, 0) = 1.;
	x_input(4, 0) = Ppiston;

#if defined(_LINUX) | (_WIN32)
	if (struc_ppts.cont_rom == false)
	{
		auto locals = py::dict("t"_a = 0, "x_input"_a = x_input, "drom"_a = drom);
		py::exec(R"(
			from rom_am import EDMD, ROM
			drom = locals()["drom"]
			pred = drom.predict(t = locals()["t"], x_input = locals()["x_input"])
		)",
				 py::globals(), locals);
		py::module_ np = py::module_::import("numpy");
		pred_ob = np.attr("real")(locals["pred"]);
	}
	else
	{
		// py::dict globals = py::globals();
		auto globals = py::dict("t"_a = 0, "x_input"_a = x_input, "erom"_a = drom, "dt"_a = Delta_t, "u_dot_t"_a = u_dot_t);
		py::exec(R"(
			from rom_am import EDMD
			import numpy as np
			from scipy.integrate import solve_ivp
			erom = locals()["erom"]
			dt = locals()["dt"]
			x_input = locals()["x_input"]

			prev_u_t = x_input[0]
			prev_u_dot_t = locals()["u_dot_t"]
			pressure = x_input[-1]

			block_m = np.block([[np.zeros((1, 5)), np.array([1])], [
				erom.model.A, np.array([0])]])
			def cont_sys(y): return block_m @ np.concatenate((y[:1, :], y[:1, :]**2, y[:1, :]**3, np.ones(
				(1, y.shape[1])), pressure * np.ones((1, y.shape[1])), y[1:2, :])).reshape((-1, y.shape[1]))
			def f(t, y): return cont_sys(y)

			sol = solve_ivp(f, [0, dt], np.array([prev_u_t, prev_u_dot_t]), t_eval=np.array([0, dt]), vectorized = True)
			pred = np.array([[sol.y[0, -1]], [sol.y[1, -1]], [cont_sys(sol.y)[-1, -1]]])
		)",
				 globals, globals);
		pred_ob = globals["pred"];
	}

	pred = pred_ob.cast<VectorXf>();
#endif
	u_t = pred(0, 0);
	u_dot_t = pred(1, 0);
	u_double_dot_t = pred(2, 0);
}

void STRUC::solve(float Delta_t)
{
	float inter;
	if (struc_ppts.rom_in_struc)
	{
		rom_model_solve(Delta_t);
	}
	else
	{
		if (struc_ppts.sdim == 1)
		{
			if (struc_ppts.spring_model == "linear")
			{
				lin_1D_model_solve(Delta_t, newm_beta, newm_gamma);
			}
			else if (struc_ppts.spring_model == "nonlinear")
			{
			}
			u_t = u_n(0);
			u_dot_t = u_dt_n(0);
			u_double_dot_t = u_ddt_n(0);
		}

		else
		{
			if (struc_ppts.spring_model == "linear")
			{
				lin_model_solve(Delta_t);
			}
			else if (struc_ppts.spring_model == "nonlinear")
			{
				nonlin_model_solve(Delta_t);
			}
			u_t += delta_u;
			inter = 4 / std::pow(Delta_t, 2) * delta_u - 4 / Delta_t * u_dot_t -
					u_double_dot_t;
			u_dot_t = u_dot_t + Delta_t / 2 * (u_double_dot_t + inter);
			u_double_dot_t = inter;
		}
	}
}

void STRUC::initialize(float presPist, Mesh mesh)
{
	// Initialisation of mesh, mass and stiffeness matrices (1D case)
	if (struc_ppts.sdim == 1)
	{
		msh = mesh;
		rigid.resize(msh.nnt, msh.nnt);
		mass.resize(msh.nnt, msh.nnt);
		assemble();

		set_BC_essential(msh.nnt - 1); // 0 Dirichlet Cond at the right boundary
		rhs = VectorXf::Zero(msh.nnt);
		rhs_term(-presPist);

		is_there_cholG = false;
		if (struc_ppts.ch_alph)
		{
			// Initialization of generalized alpha params
			ch_alpha_m = struc_ppts.ch_alpha_m;
			ch_alpha_f = struc_ppts.ch_alpha_f;
			ch_beta = struc_ppts.ch_beta;
			ch_gamma = struc_ppts.ch_gamma;
		}
		else
		{
			// Initialization of newmark params
			newm_beta = struc_ppts.newm_beta;
			newm_gamma = struc_ppts.newm_gamma;
		}
	}

	// Initialisation of the displacement
	if (struc_ppts.sdim == 0)
	{
		u_t = struc_ppts.U_0;
	}
	else
	{
		u_n = VectorXf::Zero(msh.nnt);
		for (int i = 0; i < msh.nnt; i++)
		{
		u_n(i) = struc_ppts.U_0 * std::cos(std::acos(-1) * i / (2 * (msh.nnt - 1)));
		}
		//u_n = (struc_ppts.U_0 * VectorXf::LinSpaced(Sequential, msh.nnt, 1, 0));

	}

	// Initialisation of the velocity
	if (struc_ppts.sdim == 0)
	{
		u_dot_t = 0;
	}
	else
	{
		u_dt_n = VectorXf::Zero(msh.nnt);
	}

	// Initialisation of the acceleration
	if (struc_ppts.sdim == 0)
	{
		if (struc_ppts.spring_model == "linear")
		{
			float vfg0;
			vfg0 = (presPist - 0 * struc_ppts.pres_init0) * struc_ppts.A;
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
				u_double_dot_t = (struc_ppts.A * presPist + struc_ppts.vprel[0] * (struc_ppts.u0 - u_t) + struc_ppts.mu * std::pow(struc_ppts.u0 - u_t, 3)) /
								 struc_ppts.vprel[1];
			}
		}
	}
	else
	{
		float E = struc_ppts.young, m = struc_ppts.rho_s;
		Eigen::SimplicialCholesky<SparseMatrix<float>> chol;
		chol.compute(mass);
		u_ddt_n = chol.solve(rhs - rigid * u_n);
	}

#if defined(_LINUX) | (_WIN32)
	if (struc_ppts.rom_in_struc)
	{
		py::dict globals = py::globals();

		if (struc_ppts.cont_rom == false)
		{
			py::exec(R"(
				import rom_am
				import numpy as np		
				train_data_disp = np.load("../References/discrete/train_data_disp.npy")
				train_data_veloc = np.load("../References/discrete/train_data_velocity.npy")
				train_data_accel = np.load("../References/discrete/train_data_accel.npy")
				train_data_pres = np.load("../References/discrete/train_data_pres.npy")
				train_dt = np.load("../References/discrete/train_dt.npy")
				dt = train_dt[0]
				X = train_data_disp[:1, :-1]
				Y = np.vstack((train_data_disp[:1, 1::], train_data_veloc[:1, 1::], train_data_accel[:1, 1::]))
				observables = {"X" : [lambda x : x, lambda x : x**2, lambda x : x**3, lambda x : np.ones((1, X.shape[1])), lambda x : train_data_pres[:, :-1]],  
							"Y" : [lambda x : x]}
				edmd = rom_am.EDMD()
				drom = rom_am.ROM(edmd)
				drom.decompose(X,  Y = Y, dt = dt, observables=observables)
			)",
					 globals, globals);

			py::object dt_ob = globals["dt"];
			dt_export = dt_ob.cast<float>();
		}
		else
		{
			py::exec(R"(
				import rom_am
				import numpy as np		
				train_disp = np.load("../References/continuous/train_disp.npy")
				train_accel = np.load("../References/continuous/train_accel.npy")
				train_pres = np.load("../References/continuous/train_pres.npy")
				X = train_disp[:1, :].reshape((1, -1))
				Y = train_accel
				observables = {"X" : [lambda x : x, lambda x : x**2, lambda x : x**3, lambda x : np.ones((1, X.shape[1])), lambda x : train_pres],  
								"Y" : [lambda x : x]}
				edmd = rom_am.EDMD()
				drom = rom_am.ROM(edmd)
				drom.decompose(X,  Y = Y, dt = 0.1, observables=observables)
			)",
					 globals, globals); // dt is not important here
			dt_export = 0.;
		}

		drom = globals["drom"];
	}
#endif
}

MatrixXf STRUC::mass_e(VectorXf x)
{
	float le;
	MatrixXf masse, interm(2, 2);

	le = x(1, 0) - x(0, 0);

	interm(0, 0) = 2;
	interm(0, 1) = 1;
	interm(1, 0) = 1;
	interm(1, 1) = 2;
	masse = le / 6 * interm;

	return masse;
}

MatrixXf STRUC::rigid_e(VectorXf x)
{
	float le;
	MatrixXf rigide, interm(2, 2);

	le = x(1, 0) - x(0, 0);

	interm(0, 0) = 1;
	interm(0, 1) = -1;
	interm(1, 0) = -1;
	interm(1, 1) = 1;
	rigide = 1 / le * interm;

	return rigide;
}

void STRUC::rhs_term(float p)
{
	float le;

	rhs(0) = p;
}

void STRUC::assemble()
{
	VectorXi elem_id;
	MatrixXi conec;
	VectorXf coor, coor_e;
	MatrixXf rigide, masse;
	float E, m;

	m = struc_ppts.rho_s;
	E = struc_ppts.young;

	conec = msh.get_conec();
	coor = msh.get_vcor();

	for (int ie = 0; ie < msh.nelt; ie++)
	{
		elem_id = conec(ie, seq(0, 1));
		coor_e = coor(elem_id.array(), 0);
		rigide = rigid_e(coor_e);
		masse = mass_e(coor_e);

		for (int i = 0; i < 2; i++)
		{
			for (int j = 0; j < 2; j++)
			{
				rigid.coeffRef(elem_id(i), elem_id(j)) += rigide(i, j);
				mass.coeffRef(elem_id(i), elem_id(j)) += masse(i, j);
			}
		}
	}
	rigid = E * rigid;
	mass = m * mass;
}

void STRUC::lin_1D_model_solve(float Delta_t, float beta, float gamma)
{
	SparseMatrix<float> G;
	VectorXf H, prev_u_ddt = u_ddt_n;

	if (struc_ppts.ch_alph)
	{
		if (!is_there_cholG)
		{
			G = (1 - ch_alpha_m) * mass + (1 - ch_alpha_f) * std::pow(Delta_t, 2) * ch_beta * rigid;
			chol_G.compute(G);
			is_there_cholG = true;
		}
		H = rhs - ((1 - ch_alpha_f) * std::pow(Delta_t, 2) * (.5 - ch_beta) * rigid + ch_alpha_m * mass) * u_ddt_n - ((1 + ch_alpha_f) * rigid) * u_n - (1 - ch_alpha_f) * Delta_t * rigid * u_dt_n;

		u_ddt_n = chol_G.solve(H);
		u_n = u_n + Delta_t * u_dt_n + std::pow(Delta_t, 2) * ((.5 - ch_beta) * prev_u_ddt + ch_beta * u_ddt_n);
		u_dt_n = u_dt_n + Delta_t * ((1 - ch_gamma) * prev_u_ddt + ch_gamma * u_ddt_n);
	}
	else if (struc_ppts.newm)
	{
		if (!is_there_cholG)
		{
			G = mass + beta * std::pow(Delta_t, 2) * rigid;
			chol_G.compute(G);
			is_there_cholG = true;
		}
		H = rhs - rigid * (u_n + Delta_t * u_dt_n + (.5 - beta) * std::pow(Delta_t, 2) * u_ddt_n);

		u_ddt_n = chol_G.solve(H);
		u_dt_n = u_dt_n + (1 - gamma) * Delta_t * prev_u_ddt + gamma * Delta_t * u_ddt_n;
		u_n = u_n + Delta_t * u_dt_n + (.5 - beta) * pow(Delta_t, 2) * prev_u_ddt + beta * pow(Delta_t, 2) * u_ddt_n;
	}
	else
	{
		if (!is_there_cholG)
		{
			G = mass;
			chol_G.compute(G);
			is_there_cholG = true;
		}
		H = rhs - rigid * u_n - Delta_t * rigid * u_dt_n;

		u_ddt_n = chol_G.solve(H);
		u_dt_n = u_dt_n + Delta_t * u_ddt_n;
		u_n = u_n + Delta_t * u_dt_n;
	}
}

void STRUC::store_test_data(vector<VectorXf, aligned_allocator<VectorXf>> &histo_deformation,
							vector<VectorXf, aligned_allocator<VectorXf>> &histo_dt,
							vector<VectorXf, aligned_allocator<VectorXf>> &histo_ddt)
{
	histo_deformation.push_back(u_n);
	histo_dt.push_back(u_dt_n);
	histo_ddt.push_back(u_ddt_n);
}
