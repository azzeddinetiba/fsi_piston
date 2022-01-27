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
	string spring_model;
	int nln_order;
};

inline properties load_ppts();

#endif
