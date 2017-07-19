//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   Shift/mc_tallies/FETs/FET.hh
 * \author Ryan H. Stewart
 * \date   Wed June 27 05:08:30 2017
 * \brief  FET definition.
 */
//---------------------------------------------------------------------------//


#include<iostream>
#include<vector>
#include<cstdlib>
#include<cmath>

class legendre_info
{
public:
float min;
float max;
int M;
int N;
int n_counter;
double alpha_n;
std::vector<float> ortho_const_n;
std::vector<float> ortho_const_m;
std::vector<float> A_n;
std::vector<float> A_m;
std::vector<float> a_m;
std::vector<std::vector<float> > A_n_m;
std::vector<std::vector<float> > sigma_a_n_a_m;
std::vector<float> a_hat_n;
std::vector<float> a_hat_m;
std::vector<float> a_hat_n_m;
std::vector<float> current;
std::vector<float> current_unc;
std::vector<float> var_a_n;

double scale (double x, legendre_info basis);
float rescale(float var, legendre_info basis);
void get_current (legendre_info &basis);
float get_ortho_const(int n, legendre_info & basis);
void initalize (legendre_info &basis);
};

//Initalize the info for the legendre polynomial structure
/*legendre_info::legendre_info()
{
	basis.min = -1;
	basis.max = 1;
	basis.M = 3;
	basis.N = 100;
	basis.A_n_m.resize(basis.M);
	basis.sigma_a_n_a_m.resize(basis.M);
	for(int j=0; j<basis.M; j++)
	{
		basis.A_n.push_back(0);
		basis.A_m.push_back(0);
		basis.a_hat_n.push_back(0);
		basis.a_hat_m.push_back(0);
		basis.a_hat_n_m.push_back(0);
		basis.sigma_a_n_a_m[j].resize(basis.N);
		basis.A_n_m[j].resize(basis.N);
		basis.a_m.push_back(0);
		basis.current.push_back(0);
		basis.ortho_const_n.push_back(0);
		basis.ortho_const_m.push_back(0);
		basis.current_unc.push_back(0);
		basis.var_a_n.push_back(0);
	}
}*/




class particle_info
{
public:
std::vector<float> a_n;
std::vector<float> sigma_a_n;
double b_weight;
int k_particle;
bool b_alive;
double x_tild;

void get_particle (legendre_info &basis, particle_info &a);
};

class tally_info
{


};

double Pn(int n, double x);
void basis_eval (legendre_info &basis, particle_info &a);
void get_A (legendre_info &basis, particle_info &a);

void get_a_hat (legendre_info &basis, particle_info &a);

