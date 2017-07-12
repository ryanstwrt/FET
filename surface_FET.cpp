//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   Shift/mc_tallies/FETs/surface_FET.cpp
 * \author Ryan H. Stewart
 * \date   Wed June 28 04:30:47 2017
 * \brief  Funciton for scaling to legendre space.
 */
//---------------------------------------------------------------------------//

#include<cmath>
#include"FET.hh"

using namespace std;

//Dummy random number generator for testing
double random_num ()
{
	return rand()/(double) RAND_MAX;
}

//Scales the original phase space down to Legendre phase space [-1,1]
float scale (float x, 
	     legendre_info basis)
{

	if(x < basis.min || x > basis.max)
	{
		std::cout<< "The domain does not encompass the entire range of the problem."<<endl;
	}
	else
	{
		float x_tild = 2 * ((x - basis.min)/(basis.max - basis.min))
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
void basis_eval (legendre_info &basis, 
		 particle_info &a)
{
	double x;
	for(int n = 0; n<basis.N; n++)
	{

	x=exp(- random_num());
	cout<<x<<endl;
		//calculate the legendre coefficients up to truncation 			value M for a_n and a_m
		int k=0;
		do
		{
			for(int m=0; m<basis.M; m++)
			{

				a.vec_tild[m] = scale(x, basis);
				basis.alpha_n[m] = Pn(m,a.vec_tild[m]) * 					a.b_weight;
				basis.a_m[m] = Pn(m,a.vec_tild[m]);
				a.a_n[m] += basis.alpha_n[m];
			}
			k++;
		}
		while( k < a.k_particle);
	}
}

//Increments A for each particle that passes that contribues to the current on the give surface
void get_A (legendre_info &basis, 
	    particle_info &a)
{
//Calculate A_n and A_n_m
	for(int i=0;i<basis.M;i++)
	{
		basis.A_n[i] += a.a_n[i];
		basis.A_m[i] += basis.a_m[i];
		basis.A_n_m[i] += (a.a_n[i] * basis.a_m[i]);
	}
}


//This is a dummy function until I figure out where to get the particle
void get_particle (particle_info &a)
{
	a.b_weight = random_num();
	a.k_particle = 5;
}

//Initalize the info for the legendre polynomial structure
void initalize (legendre_info &basis, 
		particle_info &a)
{
	basis.min = 0;
	basis.max = 10;
	basis.M = 5;
	basis.N = 1000;

	for(int j=0; j<basis.M; j++)
	{
		basis.A_n.push_back(0);
		basis.A_n_m.push_back(0);
		basis.A_m.push_back(0);
		basis.a_hat_n.push_back(0);
		basis.a_hat_m.push_back(0);
		basis.a_hat_n_m.push_back(0);
		basis.sigma_a_n_a_m.push_back(0);
		basis.a_m.push_back(0);
		basis.alpha_n.push_back(0);
		basis.current.push_back(0);
		basis.ortho_const_n.push_back(0);
		basis.ortho_const_m.push_back(0);
		basis.current_unc.push_back(0);
		a.a_n.push_back(0);
		a.vec_tild.push_back(0);
	}
}


int main ()
{

legendre_info basis;
particle_info a;

initalize (basis, a);
get_particle(a);


basis_eval(basis, a);
get_A (basis, a);

get_a_hat(basis, a);
basis.get_current(basis);

std::cout<<"The scaled current is described as: "<<endl;
for(int n=0; n<basis.M; n++)
{
if(n==0)
{
std::cout<<basis.current[n]<<" +/- "<<basis.current_unc[n]<<endl;
}
else
{
std::cout<<basis.current[n]<<"*x^"<<n<<" +/- "<<basis.current_unc[n]<<endl;
}
}
std::cout<<endl;
}
