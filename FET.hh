//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   Shift/mc_tallies/FETs/FET.hh
 * \author Ryan H. Stewart
 * \date   Wed June 27 05:08:30 2017
 * \brief  Definition for discrete event FETs.
 */
//---------------------------------------------------------------------------//


#include<iostream>
#include<vector>
#include<cstdlib>
#include<cmath>
//#include<chrono>

class tally_info
{
  public:
//inline tally_info ();
    int surface_index;
    int num_surfaces;
    std::vector<float> surface_indices;
    std::vector<std::vector<std::vector<float> > > surface_tallies;
    std::vector<std::vector<float> > current_matrix;
    std::vector<std::vector<float> > unc_matrix;
    std::vector<std::vector<float> > R_sqr_value;
};

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
    std::vector<float> A_n;
    std::vector<float> a_hat_n;
    std::vector<float> total_current;
    std::vector<float> current_unc;
    std::vector<float> ortho_const;
    std::vector<float> var_a_n;
    std::vector<float> unc_term1;
};

//Constructor for legendre_info
legendre_info::legendre_info ()
		 : min(0)
		  ,max(125)
		  ,M(5)
		  ,N(1000000)
		  ,n_counter(0)
{
    for(int m=0; m<M; m++)
    {
	A_n.push_back(0);
	a_hat_n.push_back(0);
	total_current.push_back(0);
	current_unc.push_back(0);
	ortho_const.push_back(0);
	var_a_n.push_back(0);
	unc_term1.push_back(0);
    }
}

//This class will either absorb the values coming out of shift, or disappear once integrated
class particle_info
{
  public:
    std::vector<float> a_n;
    double b_weight;
    int k_particle;
    bool b_alive;
    double x_tild;
    double particle_surface;

    void get_particle (legendre_info &basis, particle_info &a);
};

    void surface_eval (legendre_info &basis, particle_info &a, tally_info &tally);
    void get_A (legendre_info &basis, particle_info &a, tally_info &tally);
    void get_current (legendre_info &basis, particle_info &a, tally_info &tally);
    void get_legendre_coefficient (legendre_info &basis, particle_info &a, tally_info &tally);
    void initialize_tally_info (tally_info &tally, legendre_info &basis);
    double Pn(int n, double x);
    double scale (double x, legendre_info basis);
    float rescale(float var, legendre_info basis);
    float get_ortho_const(int n, legendre_info & basis);

