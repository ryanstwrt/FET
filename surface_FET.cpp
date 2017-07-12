//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   Shift/mc_tallies/FETs/surface_FET.cpp
 * \author Ryan H. Stewart
 * \date   Wed June 28 04:30:47 2017
 * \brief  Funciton for scaling to legendre space.
 */
//---------------------------------------------------------------------------//

#include<iostream>
#include<vector>
#include<cstdlib>
#include"FET.hh"

using namespace std;

float scale (float x, legendre_info basis)
{

	if(x < basis.min || x > basis.max)
	{
	std::cout<< "The domain does not encompass the entire range of the problem."<<endl;
	}
	else
	{
	float x_tild = 2.0 * ((x - basis.min)/(basis.max - basis.min)) - 1.0;
	return x_tild;
	}
}

float rescale (float x_tild, legendre_info basis)
{

	float x = ((x_tild+1)*(basis.max-basis.min)/2)+basis.min;;
	return x;

}

void basis_eval (legendre_info &basis, particle_info &a)
{
	double x;
	for(int n = 0; n<basis.N; n++)
	{
	x = rand() / (double) RAND_MAX;

		//calculate the legendre coefficients up to truncation 			value M for a_n and a_m
		int k=0;
		do
		{
		for(int m=0; m<basis.M; m++)
		{

		a.vec_tild[m] = scale(x, basis);

		basis.alpha_n[m] = Pn(m,a.vec_tild[m]) * 			a.b_weight;
		basis.a_m[m] = Pn(m,a.vec_tild[m]);
		a.a_n[m] += basis.alpha_n[m];
		}
		k++;
		}
		while( k < a.k_particle);
	}
}

void get_A (legendre_info &basis, particle_info &a)
{
//Calculate A_n and A_n_m
	for(int i=0;i<basis.M;i++)
	{
	basis.A_n[i] += a.a_n[i];
	basis.A_m[i] += basis.a_m[i];
	basis.A_n_m[i] += (a.a_n[i] * basis.a_m[i]);
	}
}




void get_uncertainty(legendre_info &basis)
{

}

//This is a dummy function until I figure out where to get the particle
void particle_info::get_particle (particle_info &a)
{
a.b_weight = 0.5;
a.k_particle = 5;
}

//Initalize the info for the legendre polynomial
void initalize (legendre_info &basis, particle_info &a)
{
	basis.min=0.0;
	basis.max=5.0;
	basis.M = 5;
	basis.N = 100;
	for(int j=0; j<basis.M; j++)
	{
	basis.A_n.push_back(0.0);
	basis.A_n_m.push_back(0.0);
	basis.A_m.push_back(0.0);
	basis.a_hat_n.push_back(0.0);
	basis.a_hat_m.push_back(0.0);
	basis.a_hat_n_m.push_back(0.0);
	basis.sigma_a_n_a_m.push_back(0.0);
	basis.a_m.push_back(0.0);
	basis.alpha_n.push_back(0.0);
	basis.current.push_back(0.0);
	basis.ortho_const_n.push_back(0.0);
	basis.ortho_const_m.push_back(0.0);
	basis.current_unc.push_back(0.0);
	a.a_n.push_back(0.0);
	a.vec_tild.push_back(0.0);
	}
}


int main ()
{

legendre_info basis;
particle_info a;

initalize (basis, a);
a.get_particle(a);


basis_eval(basis, a);
get_A (basis, a);

get_a_hat(basis, a);
get_current(basis);

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
