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

enum dimension
{
one,
two,
three
};

float transform (float var, legendre_info basis)
{

	if(var < basis.min || var > basis.max)
	{
	std::cout<< "The domain does not encompass the entire range of the problem."<<endl;
	}
	else
	{
	float x_tild = 2.0 * ((var - basis.min)/(basis.max - basis.min)) - 1.0;
	return x_tild;
	}
}

void basis_eval (legendre_info &basis, particle_info &a)
{
	float x;
	for(int n = 0; n<basis.N; n++)
	{
	x = rand() / (float) RAND_MAX;

		//initialize the a_n vector if not already initialized
		if(a.a_n.empty() || basis.a_m.empty() == 1)
		{
			for(int l=0; l<basis.M; l++)
			{

			}
		}
	//scale our dimension into legendre space from [-1,1] 

		//calculate the legendre coefficients up to truncation value M for a_n and a_m
		int k=0;
		do
		{
		for(int m=0; m<basis.M; m++)
		{

		a.vec_tild[m] = transform(x, basis);
		
		std::cout<<"Vec_Tild is "<<a.vec_tild[m]<<endl;
		basis.alpha_n[m] = Pn(m,a.vec_tild[m]);
		std::cout<<"alpha is "<<basis.alpha_n[m]<<endl;
		basis.alpha_n[m] *= a.weight;
		basis.a_m[m] = Pn(m,a.vec_tild[m]);
		std::cout<<"a_m is "<<basis.a_m[m]<<endl;
		a.a_n[m] += basis.alpha_n[m];
		
		std::cout<<"a_"<<n<<"_"<<m<<" = "<<a.a_n[m]<<endl;
		std::cout<<"a_"<<n<<"_"<<m<<" = "<<basis.a_m[m]<<endl;
		std::cout<<endl;
		}
		k++;
		}
		while( k < a.k_particle);
		std::cout<<"Particle "<<n<<" crossed the surface "<<k<<" times"<<endl;
		std::cout<<endl;
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


//This is a dummy function until I figure out where to get the particle
void get_particle (particle_info &a)
{
a.weight = 0.5;
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
	a.a_n.push_back(0.0);
	basis.a_m.push_back(0.0);
	basis.alpha_n.push_back(0.0);
	a.vec_tild.push_back(0.0);
	}
}

int main ()
{

legendre_info basis;
particle_info a;

//initialize the FET program
initalize (basis, a);
get_particle(a);


//Calculate the basis functions for particle a_n
basis_eval(basis, a);
//Caclulate the A_n
get_A (basis, a);

std::cout<<"This is a_n and a_m"<<endl;
for(int n = 0; n < basis.M; n++)
{
std::cout<<a.a_n[n]<<"  ";
std::cout<<basis.a_m[n]<<endl;
}
std::cout<<endl;
std::cout<<"This is A_n and A_n_m"<<endl;
for(int n = 0; n < basis.M; n++)
{
std::cout<<basis.A_n[n]<<"  ";
std::cout<<basis.A_n_m[n]<<endl;
}

get_a_hat(basis, a);

for(int n = 0; n < basis.M; n++)
{
std::cout<<basis.a_hat_n[n]<<"  ";
std::cout<<basis.a_hat_m[n]<<"  ";
std::cout<<basis.a_hat_n_m[n]<<endl;
std::cout<<endl;
std::cout<<"The "<<n<<" coefficient is "<<basis.a_hat_n[n]<<" with a std. dev. "<<basis.sigma_a_n_a_m[n]<<endl;
std::cout<<endl;
}

}
