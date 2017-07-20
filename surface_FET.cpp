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
		double x_tild = 2 * ((x - basis.min)/(basis.max - basis.min))
		 - 1;
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
	for(int n = 0; n<basis.N; n++)
	{
		x=20*random_num()-10;
		a.get_particle(basis, a);
		//calculate the legendre coefficients up to truncation value M for a_n and a_m for one particle
		//k is the number of times the particle crosses the specified surface, while m is the legendre coefficient
		for(int k=1; k<=a.k_particle; k++)
		{
			for(int m=0; m<basis.M; m++)
			{
				
				a.x_tild = scale(x, basis);
				basis.alpha_n = a.b_weight * Pn(m,a.x_tild);
				if(k==1)
				{
					
					a.a_n[m] = 0;
					a.a_n[m] = basis.alpha_n;
					basis.a_m[m] = a.a_n[m];
				}
				else
				{
					a.a_n[m] += basis.alpha_n;
					basis.a_m[m] = a.a_n[m];
				}

			}

		}

		get_A (basis, a, tally);
		basis.n_counter++;

	}
	get_a_hat(basis, a);
	get_current(basis);
}

//Increments A for each particle that passes that contribues to the current on the give surface
//verified 7/18/17
void get_A (legendre_info &basis, 
	    particle_info &a,
	    tally_info &tally)
{
//Calculate A_n and A_n_m
	for(int i=0;i<basis.M;i++)
	{
		basis.A_n[i] += a.a_n[i];
		basis.A_m[i] += basis.a_m[i];
		basis.A_n_m[i][basis.n_counter] = a.a_n[i] * basis.a_m[i];
		tally.surface_tallies[i][basis.n_counter][tally.surface_index] = basis.A_n[i];
		std::cout<<tally.surface_tallies[i][basis.n_counter][tally.surface_index]<<std::endl;
	}
}


//This is a dummy function until I figure out where to get the particle
void particle_info::get_particle (legendre_info &basis,
				  particle_info &a)
{
	a.b_weight = random_num();
	a.k_particle = 6 * random_num();
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

tally.num_surfaces = 2;
std::vector<std::vector<std::vector<float> > > surface_tallies;

tally.surface_tallies.resize(basis.M);

	for(int m=0;m<basis.M;m++)
	{
	tally.surface_tallies[m].resize(basis.N);
		for(int n=0;n<basis.N;n++)
		{
			tally.surface_tallies[m][n].resize(tally.num_surfaces);
		}
	}
}
