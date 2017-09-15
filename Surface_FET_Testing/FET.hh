
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
    inline tally_info (std::size_t poly_order);
    virtual ~tally_info() = default;

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

    std::vector<std::vector<std::vector<double> > > flux_matrix;
    std::vector<double> flux_unc;
    std::vector<double> flux_R_sqr;

    double current_total_unc;
    double flux_total_unc;


  private:
    std::size_t order;
    std::size_t terms;
};

tally_info::tally_info (std::size_t poly_order) 
			: order(poly_order)
			,terms(poly_order+1)
{

    current_matrix_x.resize(terms, 0.0);
    current_matrix_y.resize(terms, 0.0);

    flux_matrix_x.resize(terms, 0.0);
    flux_matrix_y.resize(terms, 0.0);
    flux_matrix_z.resize(terms, 0.0);

    current_unc_matrix.resize(terms, 0.0);
    current_R_sqr_value.resize(terms, 0.0);

    flux_unc_matrix.resize(terms, 0.0);
    flux_R_sqr_value.resize(terms, 0.0);

    flux_matrix.resize(terms);

   for(int m=0; m<terms; ++m)
    {
	flux_matrix[m].resize(terms);
      for(int n=0; n<terms; ++n)
      {
	flux_matrix[m][n].resize(terms);
       for(int i=0; i<terms; ++i)
       {
       flux_matrix[m][n][i] = 0;
       }
     }
    }

    flux_unc.resize(terms,0.0);
    flux_R_sqr.resize(terms,0.0);
}
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
    std::vector<double> B_n_x;
    std::vector<double> B_n_y;
    std::vector<double> B_n_z;
    std::vector<double> P_n;
    std::vector<double> b_n_x;
    std::vector<double> b_n_y;
    std::vector<double> b_n_z;
    std::vector<std::vector<std::vector<double> > > b;
    std::vector<std::vector<std::vector<double> > > B;

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
    A_n_x.resize(2*terms, 0.0);
    A_n_y.resize(2*terms, 0.0);
    B_n_x.resize(2*terms, 0.0);
    B_n_y.resize(2*terms, 0.0);
    B_n_z.resize(2*terms, 0.0);

    b_n_x.resize(terms, 0.0);
    b_n_y.resize(terms, 0.0);
    b_n_z.resize(terms, 0.0);
    b.resize(terms);
    B.resize(terms);

   for(int m=0; m<terms; ++m)
    {
    b[m].resize(terms);
    B[m].resize(terms);
      for(int n=0; n<terms; ++n)
      {
    b[m][n].resize(terms);
    B[m][n].resize(terms);
       for(int i=0; i<terms; ++i)
       {
    b[m][n][i] = 0;
    B[m][n][i] = 0;
       }
     }
    }
	x_basis[0] = -10;
	y_basis[0] = -10;
	z_basis[0] = -10;
	x_basis[1] = 10;
	y_basis[1] = 10;
	z_basis[1] = 10;
}

//This class will either absorb the values coming out of shift, or disappear once integrated
class particle_info
{
  public:
    std::vector<double> a_n;
    double b_weight;
    int k_particle;
    bool b_alive;
    double x;
    double y;
    double z;    
    double particle_surface;
    double xs_tot;
    double size;

    void get_particle (legendre_info &basis, particle_info &a);
};

class FET_solver
{
  public:

    void surface_eval (legendre_info &basis, particle_info &a, std::size_t poly_terms);
    void collision_eval (legendre_info &basis, particle_info &a, std::size_t poly_terms);
    void collision_eval2 (legendre_info &basis, particle_info &a, std::size_t poly_terms);
    void get_current (legendre_info &basis, tally_info &tally, std::size_t poly_terms, std::size_t N);
    double scale (double x, legendre_info basis);
};

