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

using namespace std;

template <class legendre_0>
//for n=0;
double P0(legendre_0 x)
{
	return 1.0;
}

//for n=1
template <class legendre_1>
legendre_1 P1(legendre_1 x)
{
	return x;
}

//for n
template <class legendre_n>
legendre_n Pn(int n, legendre_n x)
{
	if (n==0)
	{
	return P0(x);
	}

	else if (n==1)
	{
	return P1(x);
	}
//calls out the special case of x=1
	if (x == 1.0)
	{
	return 1.0;
	}

	if (x == -1.0)
	{
	return ((n % 2 ==0) ? 1.0 : -1.0);
	}
	
	if ((x==0) && (n % 2))
	{
	return 0.0;
	}

	legendre_n pn;
	pn = ((2*n-1)*x*Pn(n-1,x)-(n-1)*Pn(n-2,x))/n;

	return pn;	
}

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
		
		basis.alpha_n[m] = Pn<double>(m,a.vec_tild[m]);
		basis.alpha_n[m] *= a.b_weight;
		basis.a_m[m] = Pn<double>(m,a.vec_tild[m]);
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


//This is a dummy function until I figure out where to get the particle
void get_particle (particle_info &a)
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
	a.a_n.push_back(0.0);
	basis.a_m.push_back(0.0);
	basis.alpha_n.push_back(0.0);
	a.vec_tild.push_back(0.0);
	basis.current.push_back(0.0);
	basis.ortho_const.push_back(0.0);
	}
}

void get_current (legendre_info &basis)
{
	char x;
	for(int n=0; n<basis.M; n++)
	{
	basis.ortho_const[n] = (2.0*n+1.0)/(basis.max-basis.min);
	basis.current[n] = basis.a_hat_n[n] * basis.ortho_const[n] * Pn(n,x);
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
get_current(basis);

for(int n = 0; n < basis.M; n++)
{
std::cout<<"The "<<n<<" coefficient is "<<basis.a_hat_n[n]<<" with a std. dev. "<<basis.sigma_a_n_a_m[n]<<endl;
std::cout<<endl;
}

std::cout<<"The current is described as: "<<endl;
for(int n=0; n<basis.M; n++)
{
std::cout<<"+"<<basis.current[n]<<"*x^"<<n;
}
std::cout<<endl;

}
