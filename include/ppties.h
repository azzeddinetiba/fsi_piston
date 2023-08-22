#ifndef PPTIES_H_
#define PPTIES_H_

#include <vector>
#include <Eigen/Dense>
#include <Eigen/StdVector>

using namespace Eigen;
using namespace std;

struct properties
{
	vector<float> vprel;
	float Coeff;
	float Lspe;
	float Lsp0;
	float L_0;
	float A;
	float U_0;
	float L_t;
	float gam;
	float gamm1;
	float R;
	float C_v;
	float pres_init0;
	float rho_init0;
	float pres_init;
	float rho_init;
	float temp_init;
	float p_ext;
	float u_init;
	float e_init;
	float temp_init0;
	float umax;
	float u0;
	float mu;
	float dt;
	float young;
	float rho_s;
	bool newm;
	float newm_beta;
	float newm_gamma;
	string spring_model;
	int nln_order;
	int sdim;
	bool rom_in_struc;
	bool cont_rom;
	bool ch_alph;
	float ch_rho;
	float ch_alpha_m;
	float ch_alpha_f;
	float ch_beta;
	float ch_gamma;
	bool qs_static;
	bool amort;
	float damp;
};

inline properties load_ppts();

#endif
