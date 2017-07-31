//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   Shift/mc_tallies/FETs/surface_FET.cpp
 * \author Ryan H. Stewart
 * \date   Wed June 28 04:30:47 2017
 * \brief  Funciton for scaling to legendre space.
 */
//---------------------------------------------------------------------------//

#include"FET.hh"

using namespace std;

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
	double x_tild = 2 * ((x - basis.min)/(basis.max - basis.min)) - 1;
	return x_tild;
    }
}

//Rescales the Legendre phase space back to the original phase space
float rescale (float x_tild, 
	       legendre_info basis)
{

    float x = ((x_tild+1)*(basis.max-basis.min)/2)+basis.min;
    return x;

}

//Calculated the individual contribution from one particle, which can contribute multiple times depending on how many times it crosses the surface
//Verified 7/18/17
void surface_eval (legendre_info &basis, 
		   particle_info &a, 
		   tally_info &tally)
{

    double x;
    x=random_num();//random_num(); 

    //calculate the legendre coefficients up to truncation value M for a_n and a_m for one particle
    //k is the number of times the particle crosses the specified surface, while m is the legendre coefficient
    for(int k=1; k<=a.k_particle; k++)
    {
	for(int m=0; m<basis.M; m++)
	{
	    a.x_tild = scale(x, basis);

	    basis.alpha_n = a.b_weight * Pn(m,a.x_tild);

	    tally.surface_tallies[m][basis.n_counter][tally.surface_index] += basis.alpha_n;
	}
    }
   basis.n_counter++;
}


//This is a dummy function until I figure out where to get the particle
void particle_info::get_particle (legendre_info &basis,
				  particle_info &a)
{
    a.b_weight = random_num();
    a.k_particle = 1;//4 * random_num() + 1;
    a.particle_surface = .25;//random_num();
    if(basis.n_counter==0)
    {
	for(int m=0;m<basis.M;m++)
	{
	    a.a_n.push_back(0);
	}
    }
}

//Initialize the surface tallies matrix and the surface index matrix
void initialize_tally_info (tally_info &tally, legendre_info &basis)
{

    tally.num_surfaces = 1;
    std::vector<std::vector<std::vector<float> > > surface_tallies;
    std::vector<std::vector<float> > current_matrix;
    std::vector<std::vector<float> > coefficient_matrix;
    std::vector<std::vector<float> > unc_matrix;
    std::vector<std::vector<float> > R_sqr_value;

    tally.surface_tallies.resize(basis.M);
    tally.current_matrix.resize(basis.M);
    tally.coefficient_matrix.resize(basis.M);
    tally.unc_matrix.resize(basis.M);
    tally.R_sqr_value.resize(basis.M);

    for(int m=0;m<basis.M;m++)
    {
 	tally.surface_tallies[m].resize(basis.N);
	tally.current_matrix[m].resize(tally.num_surfaces);
	tally.coefficient_matrix[m].resize(tally.num_surfaces);
	tally.unc_matrix[m].resize(tally.num_surfaces);
	tally.R_sqr_value[m].resize(tally.num_surfaces);

	for(int n=0;n<basis.N;n++)
	{
	     tally.surface_tallies[m][n].resize(tally.num_surfaces);
	}
    }
}
