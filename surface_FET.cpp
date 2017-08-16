//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   Shift/mc_tallies/FETs/surface_FET.cpp
 * \author Ryan H. Stewart
 * \date   Wed June 28 04:30:47 2017
 * \brief  Funciton for scaling to legendre space.
 */
//---------------------------------------------------------------------------//

#include"FET.hh"
#include "Distribution.hh"

using namespace std;

Distribution fluxshape;

//Dummy random number generator for testing
//Verified 7/18/17
double random_num ()
{
    return rand()/(double) RAND_MAX;
}

//Scales the original phase space down to Legendre phase space [-1,1]
//Verified 7/18/17
double scale (double x, 
	     legendre_info basis)
{
    if(x < basis.x_basis[basis.surface_index] || x > basis.x_basis[basis.surface_index+1])
    {
	std::cout<< "The domain does not encompass the entire range of the problem."<<endl;
    }
    else
    {
	return (2 * ((x - basis.x_basis[basis.surface_index])/(basis.x_basis[basis.surface_index+1] - basis.x_basis[basis.surface_index])) - 1);

    }
}

//Calculated the individual contribution from one particle, which can contribute multiple times depending on how many times it crosses the surface
//Verified 7/18/17
void surface_eval (legendre_info &basis, 
		   particle_info &a, std::size_t poly_terms)
{
    double temp_var;
    double x_tild = scale(a.b_weight, basis);
    double y_tild = scale(a.b_weight, basis);
    double ratio = fluxshape(x_tild)/ a.k_particle;
    double ratio2 = fluxshape(y_tild)/ a.k_particle;
    std::vector<double> a_n_x(poly_terms, 0.0);
    std::vector<double> a_n_y(poly_terms, 0.0);
    std::vector<double> P_n_x = basis.Pn(poly_terms, x_tild);
    std::vector<double> P_n_y = basis.Pn(poly_terms, y_tild);

    //calculate the legendre coefficients up to truncation value M for a_n and a_m for one particle
    //k is the number of times the particle crosses the specified surface, while m is the legendre coefficient
    for(int k=1; k<=a.k_particle; ++k)
    {
    	    for(int m=0; m<poly_terms; ++m)
    	    {
        	a_n_x[m] += P_n_x[m];
		a_n_y[m] += P_n_y[m];
	    }

    }
    for(int m=0; m<poly_terms; ++m)
    {
	temp_var = a_n_x[m] * ratio;
	basis.A_n[m] += temp_var;
	temp_var = a_n_y[m] * ratio2;
	basis.A_n_y[m] += temp_var;
	basis.A_m[m] += pow(temp_var,2);
    }

   basis.n_counter[0]++;
}


//This is a dummy function until I figure out where to get the particle
void particle_info::get_particle (legendre_info &basis,
				  particle_info &a)
{
    a.b_weight = (basis.x_basis[basis.surface_index+1] - basis.x_basis[basis.surface_index]) * random_num();
    a.k_particle = 1;//4 * random_num() + 1;
    a.particle_surface = .25;//random_num();
}

//Initialize the surface tallies matrix and the surface index matrix
void initialize_tally_info (tally_info &tally, std::size_t poly_order)
{
    tally.current_matrix.resize(poly_order, 0.0);
    tally.current_matrix_y.resize(poly_order, 0.0);
    tally.unc_matrix.resize(poly_order, 0.0);
    tally.R_sqr_value.resize(poly_order, 0.0);
}
