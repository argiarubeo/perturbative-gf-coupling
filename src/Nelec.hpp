//
//  Nelec.hpp
//
//
//  Created by Argia Rubeo on 02/05/17.
//
//

#ifndef Nelec_h
#define Nelec_h

#include <stdio.h>

#include <map>
#include <iomanip>
#include <cmath>
#include <iostream>
#include <fstream>
#include <armadillo>

#endif /* Nelec_h */

using namespace std;

//global variables
char DELIMITER = ' ';

int L = 0;
int T = 0;
int Tmax=72;

double x_0_T = 0.0;

//momenta
double p_1 = 0.0;
double p_2 = 0.0;
double p_3 = 0.0;
double p_0 = 0.0;

//momenta on the lattice (p hat)
double p_hat_0 = 0.0;
double p_hat_1 = 0.0;
double p_hat_2 = 0.0;
double p_hat_3 = 0.0;

// p hat squared
double p_hat_2nd = 0;
double p_hat_2nds = 0;

//definition to make shorter expression for zeuthen flow
double prefact0= 0;
double prefact1= 0;
double prefact2= 0;
double prefact3= 0;

// p hat to the fourth
double p_hat_4th = 0.0;
double p_hat_4ths = 0.0;

//product of momenta on the lattice

double p0_2nd = 0.0;
double p1_2nd = 0.0;
double p2_2nd = 0.0;
double p3_2nd = 0.0;
double p0_4th = 0.0;
double p1_4th = 0.0;
double p2_4th = 0.0;
double p3_4th = 0.0;

double p0p1 = 0.0;
double p1p0 = 0.0;
double p0p2 = 0.0;
double p2p0 = 0.0;
double p0p3 = 0.0;
double p3p0 = 0.0;
double p1p2 = 0.0;
double p2p1 = 0.0;
double p1p3 = 0.0;
double p3p1 = 0.0;
double p2p3 = 0.0;
double p3p2 = 0.0;

//gauge parameter for heat kernel (should be positive)
double alpha = 1.0;
//gauge parameter for action
double lambda = 1.0;


double Nel[200];
double firstder_Nel[200];

double El[200];
double Elpl[200];
double Elpls[200];
double Elcl[200];
double secder_Elpl[200];
double secder_Elpls[200];
double secder_Elcl[200];
double secder_Elplcl[200];
double secder_Elimp[200];
double secder_El[200];
double secder_Nel[200];
double Elplcl[200];
double Elimp[200];

double firstder_Elpls[200];
double firstder_Elcl[200];
double firstder_Elpl[200];
double firstder_Elplcl[200];
double firstder_Elimp[200];
double firstder_El[200];

double Nel_p_eps[200];
double El_p_eps[200];
double Elpl_p_eps[200];
double Elpls_p_eps[200];
double Elcl_p_eps[200];
double Elplcl_p_eps[200];
double Elimp_p_eps[200];

double Nel_m_eps[200];
double El_m_eps[200];
double Elpl_m_eps[200];
double Elpls_m_eps[200];
double Elcl_m_eps[200];
double Elplcl_m_eps[200];
double Elimp_m_eps[200];

double de_Nel[200];
double de_El[200];
double de_Elpl[200];
double de_Elpls[200];
double de_Elcl[200];
double de_Elplcl[200];
double de_Elimp[200];

//double array_cont_mag[Tmax];
//double array_cont_el[Tmax];
map <double, double> cont_el_map;


//new variables nedeed for c_b

double c_b = 0.0;
double eps = 0.0;

double c_hat_1 = 0.0;
double c_hat_2 = 0.0;
double c_hat_3 = 0.0;
double c_hat_0 = 0.0;
double c_hat_squared = 0;

double p_dot_1 = 0.0;
double p_dot_2 = 0.0;
double p_dot_3 = 0.0;
double p_dot_0 = 0.0;
double p_dot_squared = 0.0;

//continuum values

double E_cont_el_T2 = 0.0;
double E_cont_el_T4 = 0.0;

double secder_cont_el_T2 = 0.0;
