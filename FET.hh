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
inline legendre_info ();

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


//void initalize (legendre_info &basis);
};

//Initalize the info for the legendre polynomial structure
legendre_info::legendre_info ()
{
	min = -1;
	max = 1;
	M = 3;
	N = 100;
	A_n_m.resize(M);
	sigma_a_n_a_m.resize(M);
	for(int j=0; j<M; j++)
	{
		A_n.push_back(0);
		A_m.push_back(0);
		a_hat_n.push_back(0);
		a_hat_m.push_back(0);
		a_hat_n_m.push_back(0);
		sigma_a_n_a_m[j].resize(N);
		A_n_m[j].resize(N);
		a_m.push_back(0);
		current.push_back(0);
		ortho_const_n.push_back(0);
		ortho_const_m.push_back(0);
		current_unc.push_back(0);
		var_a_n.push_back(0);
	}
}




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
double scale (double x, legendre_info basis);
float rescale(float var, legendre_info basis);
void get_current (legendre_info &basis);
float get_ortho_const(int n, legendre_info & basis);
void get_a_hat (legendre_info &basis, particle_info &a);

