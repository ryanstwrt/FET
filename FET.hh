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
    inline legendre_info (std::size_t poly_order, std::size_t num_surfaces);
    virtual ~legendre_info() = default;

    double min;
    double max;
    std::vector<double> n_counter;
    std::vector<double> A_n;
    std::vector<double> A_m;
    std::vector<double> Pn(std::size_t poly_terms, double x);
    std::vector<double> P_n;

  private:
    std::size_t order;
    std::size_t terms;
    std::size_t num_surfaces;
    double load(std::size_t index) const;
    void save(std::size_t index, double value);
};

//Constructor for legendre_info
legendre_info::legendre_info (std::size_t poly_order, std::size_t num_surfaces)
		 : order(poly_order)
		  ,terms(poly_order+1)
		  ,min(0)
		  ,max(2)
		  ,n_counter(num_surfaces, 0.0)
		  ,P_n(poly_order+1, 0.0)
{
    A_n.resize(terms, 0.0);
    A_m.resize(terms, 0.0);
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

    void surface_eval (legendre_info &basis, particle_info &a, std::size_t poly_terms);
    void get_current (legendre_info &basis, tally_info &tally, std::size_t poly_terms, std::size_t N);
    void initialize_tally_info (tally_info &tally, std::size_t poly_terms);

    double scale (double x, legendre_info basis);
    double rescale (double x_tild, legendre_info basis);

