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

void STRUC::rom_model_solve(float Delta_t)
{
	VectorXf pred(3), x_input(5);
	py::object pred_ob;

	x_input(0, 0) = u_t;
	x_input(1, 0) = pow(u_t, 2);
	x_input(2, 0) = pow(u_t, 3);
	x_input(3, 0) = 1.;
	x_input(4, 0) = Ppiston;

	if (struc_ppts.cont_rom == false)
	{
		auto locals = py::dict("t"_a = 0, "x_input"_a=x_input, "drom"_a=drom);
		py::exec(R"(
			from rom_am import EDMD, ROM
			drom = locals()["drom"]
			pred = drom.predict(t = locals()["t"], x_input = locals()["x_input"])
		)", py::globals(), locals);
		py::module_ np = py::module_::import("numpy");
		pred_ob = np.attr("real")(locals["pred"]);
	}
	else
	{
		//py::dict globals = py::globals();
		auto globals = py::dict("t"_a = 0, "x_input"_a=x_input, "erom"_a=drom, "dt"_a=Delta_t, "u_dot_t"_a=u_dot_t);
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
		)", globals, globals);
		pred_ob = globals["pred"];
	}


	pred = pred_ob.cast<VectorXf>();
	u_t = pred(0, 0);
	u_dot_t = pred(1, 0);
	u_double_dot_t = pred(2, 0);
}

void STRUC::solve(float Delta_t)
{
	float inter;
	if(struc_ppts.rom_in_struc)
	{
		rom_model_solve(Delta_t);
	}
	else
	{
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
			)", globals, globals);

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
			)", globals, globals); // dt is not important here
			dt_export = 0.;
		}

		drom = globals["drom"];
	}
}
