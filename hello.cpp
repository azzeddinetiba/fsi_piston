/*
 * hello.cpp
 *
 *  Created on: Dec 10, 2021
 *      Author: tiba
 */

// Your First C++ Program

#include <iostream>
#include <Eigen/Dense>

using namespace Eigen;
using namespace std;

VectorXf temperature (float C_v, MatrixXf vsol) {


	  VectorXf temp(vsol.rows());

	  ArrayXf col0, col1, col2, res;
	  col0 = vsol(all, 0).array();
	  col1 = vsol(all, 1).array();
	  col2 = vsol(all, 2).array();

	  res = (col2/col0 - (.5 * col1.square())/col0.square()) / C_v;

	  temp = res.matrix();

	  return temp;


}


VectorXf pressure (float R, VectorXf vtemp, MatrixXf vsol) {


	  VectorXf pres(vsol.rows());

	  ArrayXf col0, res;
	  col0 = vsol(all, 0).array();


	  res = col0 * R * vtemp.array();

	  pres = res.matrix();

	  return pres;


}

void data_init() {

	int nelt, nnt, nnel, ndln, ndle, ndlt;
	float L_0, L_t, U_0, A;

	nelt = 70 ;
	nnt  = nelt+1;
	nnel = 2;
	ndln = 3;
	ndle = ndln*nnel;
	ndlt = nnt*ndln;

	L_0 = 1;
	A = 1;

	U_0 = .2;
	L_t = L_0 + U_0;

    // Vector of coordinates
	VectorXf x, vcor, vcor0, vcor_n, vcor_np1;
	vcor = VectorXf::LinSpaced(Sequential, L_t, 0, nnt);
	vcor0 = vcor;
	vcor_n = vcor;
	vcor_np1 = vcor;

	// Connectivity
	MatrixXd connec(nelt, 2);
	connec.col(0) = VectorXd::LinSpaced(Sequential, nelt, 0, nnt-2);
	connec.col(1) = VectorXd::LinSpaced(Sequential, nelt, 1, nnt-1);

	// Physical properties
	float gam, gamm1, R, C_v, pres_init0, temp_init0, rho_init0, pres_init, rho_init, temp_init, p_ext, u_init, e_init;

	gam   = 1.4;     // the specific heat ratio of the gas
	gamm1 = gam-1.;
	R     = 287;     // the individual gas constant
	C_v   = R/gamm1; // the specific heat capacity of the gas

	pres_init0 = 1E5;  // initial pressure for chamber length = L0
	temp_init0 = 300;  // initial temperature
	rho_init0  = pres_init0/gamm1/C_v/temp_init0; // initial volumic mass

    pres_init = pres_init0 * pow((L_0/L_t), gam);

    rho_init  = rho_init0 * pow((pres_init/pres_init0), (1./gam));
    temp_init = pres_init/rho_init/gamm1/C_v;
    p_ext     = 0*pres_init0; // pressure on the right of the piston

    // we set the initial fluid velocity and the initial total fluid energy
    u_init = 0.;
    e_init = pres_init/gamm1/rho_init + 0.5 * pow(u_init, 2.);

    // initialisation of the solution matrix
    MatrixXf vsol(nnt, 3);
    ArrayXXf vsol_i(nnt, 3);
    vsol_i.col(0) = rho_init * VectorXf::Ones(nnt).array();
    vsol_i.col(1) = rho_init * u_init * VectorXf::Ones(nnt).array();
    vsol_i.col(2) = rho_init * e_init * VectorXf::Ones(nnt).array();

    vsol = vsol_i.matrix();


    VectorXf vpres, vtemp, vcelerity;
    ArrayXf vcelerity_i;
    float cel_init;
    vtemp = temperature(C_v, vsol);
    vpres = pressure(R, vtemp, vsol);
    vcelerity_i = (gam * gamm1 * C_v * vtemp.array()).sqrt();
    vcelerity = vcelerity_i.matrix();
    cel_init = vcelerity(0);



    // Initialisation
    MatrixXf vflux(nnt, ndln), wx(nnt, 1), presL2t;

}


int main() {

    cout << "Hello World!";

    VectorXf temp, pres;
    MatrixXf m(2,3);
    m(0,0) = 3;
    m(1,0) = 2.5;
    m(0,1) = -1;

    temp = temperature(1., m);

    pres = pressure(1, temp, m);

    data_init();


	VectorXd x1, x2;
	x1 = VectorXd::LinSpaced(Sequential, 71-1, 1, 70);

    return 0;
}


