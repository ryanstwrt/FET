//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   Shift/mc_tallies/FETs/FET.hh
 * \author Ryan H. Stewart
 * \date   Wed June 27 05:08:30 2017
 * \brief  Definition for discrete event FETs.
 */
//---------------------------------------------------------------------------//


#include<iostream>
#include<fstream>
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
    std::vector<double> surface_indices;
    std::vector<double> current_matrix;
    std::vector<double> unc_matrix;
    std::vector<double> R_sqr_value;
};

class legendre_info
{
  public:
    inline legendre_info ();
    virtual ~legendre_info() = default;

    double min;
    double max;
    int M;
    int N;
    int n_counter;
    std::vector<double> A_n;
    std::vector<double> A_m;
    std::vector<double> a_n;
    std::vector<double> current_unc;
    std::vector<double> ortho_const;
    std::vector<double> var_a_n;
};

//Constructor for legendre_info
legendre_info::legendre_info ()
		 : min(0)
		  ,max(2)
		  ,M(6)
		  ,N(1e8)
		  ,n_counter(0)
{
    A_n.resize(M, 0.0);
    A_m.resize(M, 0.0);
    a_n.resize(M, 0.0);
    var_a_n.resize(M, 0.0);
    current_unc.resize(M, 0.0);
    ortho_const.resize(M, 0.0);
}

//This class will either absorb the values coming out of shift, or disappear once integrated
class particle_info
{
  public:
    std::vector<double> a_n;
    double b_weight;
    int k_particle;
    bool b_alive;
    double x_tild;
    double particle_surface;

    void get_particle (legendre_info &basis, particle_info &a);
};

    void surface_eval (legendre_info &basis, particle_info &a);
    void get_current (legendre_info &basis, particle_info &a, tally_info &tally);
    void initialize_tally_info (tally_info &tally, std::size_t poly_order);
    double Pn(int n, double x);
    double scale (double x, legendre_info basis);
    double rescale (double x_tild, legendre_info basis);

