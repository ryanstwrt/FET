//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   Shift/mc_tallies/FETs/FET.hh
 * \author Ryan H. Stewart
 * \date   Wed June 27 05:08:30 2017
 * \brief  FET definition.
 */
//---------------------------------------------------------------------------//


#include <iostream>
#include <vector>

struct legendre_info
{
float min;
float max;
int M;
int N;
std::vector<float> ortho_const_n;
std::vector<float> ortho_const_m;
std::vector<float> alpha_n;
std::vector<float> A_n;
std::vector<float> A_m;
std::vector<float> a_m;
std::vector<float> A_n_m;
std::vector<float> sigma_a_n_a_m;
std::vector<float> a_hat_n;
std::vector<float> a_hat_m;
std::vector<float> a_hat_n_m;
std::vector<float> current;
std::vector<float> current_unc;
};

class particle_info
{
public:
std::vector<float> a_n;
std::vector<float> vec_tild;
std::vector<float> sigma_a_n;
double b_weight;
int k_particle;
bool b_alive;
void get_particle (particle_info &a);
};

void get_particle (particle_info &a);
double Pn(int n, double x);
void basis_eval (legendre_info &basis, particle_info &a);
float scale (float var, legendre_info basis);
void get_A (legendre_info &basis, particle_info &a);
void initialize (legendre_info &basis, particle_info &a);
void get_a_hat (legendre_info &basis, particle_info &a);
void get_current (legendre_info &basis);
float get_ortho_const(int n, legendre_info & basis);


