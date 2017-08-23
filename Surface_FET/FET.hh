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
    std::vector<double> surface_indices;
    std::vector<double> current_matrix_x;
    std::vector<double> current_matrix_y;
    std::vector<double> flux_matrix_x;
    std::vector<double> flux_matrix_y;
    std::vector<double> flux_matrix_z;

    std::vector<double> current_unc_matrix;
    std::vector<double> flux_unc_matrix;
    std::vector<double> current_R_sqr_value;
    std::vector<double> flux_R_sqr_value;

    double current_total_unc;
    double flux_total_unc;
};
/*
class collision_tally_info
{
  private:
    double total_xc;

};*/

class legendre_info
{
  public:
    inline legendre_info (std::size_t poly_order, std::size_t num_surfaces);
    virtual ~legendre_info() = default;

    int surface_index;
    std::size_t num_surfaces;
    std::vector<double> x_basis;
    std::vector<double> y_basis;
    std::vector<double> z_basis;
    std::vector<double> n_counter;
    std::vector<double> A_n_x;
    std::vector<double> A_n_y;
    std::vector<double> A_m_x;
    std::vector<double> A_m_y;
    std::vector<double> B_n_x;
    std::vector<double> B_n_y;
    std::vector<double> B_n_z;
    std::vector<double> B_m_x;
    std::vector<double> B_m_y;
    std::vector<double> B_m_z;
    std::vector<double> P_n;

    std::vector<double> Pn(std::size_t poly_terms, double x);

  private:
    std::size_t order;
    std::size_t terms;
    double load(std::size_t index) const;
    void save(std::size_t index, double value);
};

//Constructor for legendre_info
legendre_info::legendre_info (std::size_t poly_order, std::size_t num_surfaces)
		 : order(poly_order)
		  ,terms(poly_order+1)
		  ,num_surfaces(num_surfaces)
		  ,n_counter(num_surfaces, 0.0)
		  ,P_n(poly_order+1, 0.0)
{
    x_basis.resize(2*num_surfaces, 0.0);
    y_basis.resize(2*num_surfaces, 0.0);
    z_basis.resize(2*num_surfaces, 0.0);
    A_n_x.resize(terms, 0.0);
    A_n_y.resize(terms, 0.0);
    A_m_x.resize(terms, 0.0);
    A_m_y.resize(terms, 0.0);
    B_n_x.resize(terms, 0.0);
    B_n_y.resize(terms, 0.0);
    B_n_z.resize(terms, 0.0);
    B_m_x.resize(terms, 0.0);
    B_m_y.resize(terms, 0.0);
    B_m_z.resize(terms, 0.0);
    for(int n=0; n<2*num_surfaces; ++n)
    {
	x_basis[n] = 2*n;
	y_basis[n] = 2*n;
	z_basis[n] = 2*n;
    }
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
    void collision_eval (legendre_info &basis, particle_info &a, std::size_t poly_terms);
    void get_current (legendre_info &basis, tally_info &tally, std::size_t poly_terms, std::size_t N);
    void initialize_tally_info (tally_info &tally, std::size_t poly_terms);
    double scale (double x, legendre_info basis);
