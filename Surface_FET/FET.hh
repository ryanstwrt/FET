
//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   Shift/mc_tallies/FETs/FET.hh
 * \author Ryan H. Stewart
 * \date   Wed June 27 05:08:30 2017
 * \mod    Thur Sept 20 04:53:00 2017
 * \brief  Definition for discrete event FETs.
 */
//---------------------------------------------------------------------------//


#include<iostream>
#include<fstream>
#include<vector>
#include<cstdlib>
#include<cmath>

//Holds the information regarding the final current and flux
class tally_info
{
  public:
    inline tally_info (std::size_t poly_order);
    virtual ~tally_info() = default;

    double R_great_10;
    double R_great_1;
    double total_coeff;
    std::vector<std::vector<double> > current_matrix;
    std::vector<std::vector<double> > current_unc_matrix;
    std::vector<std::vector<double> > current_R_matrix;

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

    current_matrix.resize(terms);
    current_unc_matrix.resize(terms);
    current_R_matrix.resize(terms);

    flux_matrix.resize(terms);
    flux_unc_matrix.resize(terms);
    flux_R_matrix.resize(terms);

//This loop initializes both the flux matrix, and the associated uncertainty matrix. The sum of the entire matrix will yield the flux for the system.
    for(int m=0; m<terms; ++m)
    {
      current_matrix[m].resize(terms);
      current_unc_matrix[m].resize(terms);
      current_R_matrix[m].resize(terms);

      flux_matrix[m].resize(terms);
      flux_unc_matrix[m].resize(terms);
      flux_R_matrix[m].resize(terms);

      for(int n=0; n<terms; ++n)
      {
        current_matrix[m][n] = 0;
        current_unc_matrix[m][n] = 0;
        current_R_matrix[m][n] = 0;

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

//This class holds all of the inforomation regarding the estimate for the coefficients
class legendre_info
{
  public:
    inline legendre_info (std::size_t poly_order, std::size_t num_tallies);
    virtual ~legendre_info() = default;

    std::size_t num_tallies;
    std::vector<double> x_basis;
    std::vector<double> y_basis;
    std::vector<double> z_basis;
    std::vector<double> n_counter;
    std::vector<double> P_n;

    std::vector<std::vector<double> >  a;
    std::vector<std::vector<double> >  A;
    std::vector<std::vector<double> >  A_unc;

    std::vector<std::vector<std::vector<double> > > b;
    std::vector<std::vector<std::vector<double> > > B;
    std::vector<std::vector<std::vector<double> > > B_unc;

    std::vector<double> Pn(std::size_t poly_terms, double x);

  private:
    std::size_t order;
    std::size_t terms;
    double load(std::size_t index) const;
    void save(std::size_t index, double value);
};

//Constructor for legendre_info
legendre_info::legendre_info (std::size_t poly_order, std::size_t num_tallies)
		 : order(poly_order)
		  ,terms(poly_order+1)
		  ,num_tallies(num_tallies)
		  ,n_counter(num_tallies, 0.0)
		  ,P_n(poly_order+1, 0.0)
{
    x_basis.resize(2*num_tallies, 0.0);
    y_basis.resize(2*num_tallies, 0.0);
    z_basis.resize(2*num_tallies, 0.0);

    a.resize(terms);
    A.resize(terms);
    A_unc.resize(terms);
    b.resize(terms);
    B.resize(terms);
    B_unc.resize(terms);

//This loop initilizes the coefficient matrix for both the current and the flux
//To Do: Create some statement to differentiate the two types of FETs (surface vs collision) and only intialize the coefficient matrix needed
   for(int m=0; m<terms; ++m)
    {
      a[m].resize(terms);
      A[m].resize(terms);
      A_unc[m].resize(terms);

      b[m].resize(terms);
      B[m].resize(terms);
      B_unc[m].resize(terms);
      for(int n=0; n<terms; ++n)
      {
        a[m][n] = 0;
        A[m][n] = 0;
        A_unc[m][n] = 0;;

        b[m][n].resize(terms);
        B[m][n].resize(terms);
        B_unc[m][n].resize(terms);
        for(int i=0; i<terms; ++i)
        {
          b[m][n][i] = 0;
          B[m][n][i] = 0;
          B_unc[m][n][i] = 0;
        }
      }
    }
//Currently sets the boundries for the x,y, and z dimensions
//To Do: Figure out a better way to store this information, perhaps a 3 x 3 x 3 matrix to be initialized above?
	x_basis[0] = -1;
	y_basis[0] = -1;
	z_basis[0] = -5;
	x_basis[1] = 1;
	y_basis[1] = 1;
	z_basis[1] = 5;
}

//This class will either absorb the values coming out of shift, or disappear once integrated
//To Do: position will be given with a number space number space number (ex 1 0 1). find a way to seperate this out
class particle_info
{
  public:
    double wt;
    bool alive;
    bool domain;
    double x;
    double y;
    double z;    
    double xs_tot;
    double size;
};

//Solver class simply groups all of the various functions needed to solve for the FET
class FET_solver
{
  public:
    void surface_eval (legendre_info &basis, particle_info &a, std::size_t poly_terms);
    void surface_eval2 (legendre_info &basis, particle_info &a, std::size_t poly_terms);
    void collision_eval (legendre_info &basis, particle_info &a, std::size_t poly_terms);
    void collision_eval2 (legendre_info &basis, particle_info &a, std::size_t poly_terms);
    void get_current (legendre_info &basis, tally_info &tally, std::size_t poly_terms);
    void cleanup (tally_info &tally, std::size_t poly_terms);
    double scale (double x, std::vector<double> dimension_basis);
};

