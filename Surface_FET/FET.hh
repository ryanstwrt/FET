
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

    std::vector<double> current_unc_matrix;
    std::vector<double> current_R_sqr_value;

    double current_total_unc;
    double flux_total_unc;

    std::vector<std::vector<std::vector<double> > > flux_matrix;
    std::vector<std::vector<std::vector<double> > > flux_unc_matrix;
    std::vector<std::vector<std::vector<double> > > flux_R_matrix;

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

    current_unc_matrix.resize(terms, 0.0);
    current_R_sqr_value.resize(terms, 0.0);

    flux_matrix.resize(terms);
    flux_unc_matrix.resize(terms);
    flux_R_matrix.resize(terms);

//This loop initializes both the flux matrix, and the associated uncertainty matrix. The sum of the entire matrix will yield the flux for the system.
   for(int m=0; m<terms; ++m)
    {
      flux_matrix[m].resize(terms);
      flux_unc_matrix[m].resize(terms);
      flux_R_matrix[m].resize(terms);

      for(int n=0; n<terms; ++n)
      {
	flux_matrix[m][n].resize(terms);
        flux_unc_matrix[m][n].resize(terms);
        flux_R_matrix[m][n].resize(terms);

        for(int i=0; i<terms; ++i)
        {
          flux_matrix[m][n][i] = 0;
          flux_unc_matrix[m][n][i] = 0;
          flux_R_matrix[m][n][i] = 0;
        }
      }
    }


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
    std::vector<double> P_n;
    std::vector<std::vector<std::vector<double> > > a;
    std::vector<std::vector<std::vector<double> > > A;
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
    a.resize(terms);
    A.resize(terms);


    b.resize(terms);
    B.resize(terms);

//This loop initilizes the coefficient matrix for both the current and the flux
//To Do: Create some statement to differentiate the two and only intialize the coefficient matrix needed
   for(int m=0; m<terms; ++m)
    {
      a[m].resize(terms);
      A[m].resize(terms);
      b[m].resize(terms);
      B[m].resize(terms);
      for(int n=0; n<terms; ++n)
      {
        a[m][n].resize(terms);
        A[m][n].resize(terms);
        b[m][n].resize(terms);
        B[m][n].resize(terms);
        for(int i=0; i<terms; ++i)
        {
          a[m][n][i] = 0;
          A[m][n][i] = 0;
          b[m][n][i] = 0;
          B[m][n][i] = 0;
        }
      }
    }
//Currently sets the boundries for the x,y, and z dimensions
//To Do: Figure out a better way to store this information, perhaps a 3 x 3 x 3 matrix to be initialized above?
	x_basis[0] = -1;
	y_basis[0] = -1;
	z_basis[0] = -1;
	x_basis[1] = 1;
	y_basis[1] = 1;
	z_basis[1] = 1;
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

//Solver class simply groups all of the various functions needed to solve for the FET
class FET_solver
{
  public:

    void surface_eval (legendre_info &basis, particle_info &a, std::size_t poly_terms);
    void collision_eval (legendre_info &basis, particle_info &a, std::size_t poly_terms);
    void collision_eval2 (legendre_info &basis, particle_info &a, std::size_t poly_terms);
    void get_current (legendre_info &basis, tally_info &tally, std::size_t poly_terms, std::size_t N);
    double scale (double x, legendre_info basis);
};

