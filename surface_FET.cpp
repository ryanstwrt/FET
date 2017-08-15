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
    if(x < basis.min || x > basis.max)
    {
	std::cout<< "The domain does not encompass the entire range of the problem."<<endl;
    }
    else
    {
	return (2 * ((x - basis.min)/(basis.max - basis.min)) - 1);

    }
}

//Calculated the individual contribution from one particle, which can contribute multiple times depending on how many times it crosses the surface
//Verified 7/18/17
void surface_eval (legendre_info &basis, 
		   particle_info &a, std::size_t poly_terms)
{
    double temp_var;
    double x_tild = scale(a.b_weight, basis);
    double ratio = fluxshape(x_tild) / a.k_particle;
    std::vector<double> alpha_n(poly_terms, 0.0);
    std::vector<double> a_n(poly_terms, 0.0);
    basis.Pn(poly_terms, x_tild);

    //calculate the legendre coefficients up to truncation value M for a_n and a_m for one particle
    //k is the number of times the particle crosses the specified surface, while m is the legendre coefficient
    for(int k=1; k<=a.k_particle; ++k)
    {
    	for(int m=0; m<poly_terms; ++m)
    	{
        a_n[m] += basis.P_n[m];
	}

    }
    for(int m=0; m<poly_terms; ++m)
    {
	temp_var = a_n[m] * ratio;
	basis.A_n[m] += temp_var;
	basis.A_m[m] += pow(temp_var,2);
    }

   basis.n_counter[0]++;
}


//This is a dummy function until I figure out where to get the particle
void particle_info::get_particle (legendre_info &basis,
				  particle_info &a)
{
    a.b_weight = (basis.max - basis.min) * random_num();
    a.k_particle = 1;//4 * random_num() + 1;
    a.particle_surface = .25;//random_num();
}

//Initialize the surface tallies matrix and the surface index matrix
void initialize_tally_info (tally_info &tally, std::size_t poly_order)
{
    tally.num_surfaces = 1;
    std::vector<double> current_matrix;
    std::vector<double> unc_matrix;
    std::vector<double> R_sqr_value;

    tally.current_matrix.resize(poly_order, 0.0);
    tally.unc_matrix.resize(poly_order, 0.0);
    tally.R_sqr_value.resize(poly_order, 0.0);
}
