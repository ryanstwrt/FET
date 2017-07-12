//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   Shift/mc_tallies/FETs/surface_FET.cpp
 * \author Ryan H. Stewart
 * \date   Wed June 28 04:30:47 2017
 * \brief  Funciton for scaling to legendre space.
 */
//---------------------------------------------------------------------------//

#include"FET.hh"

//Solves the the coefficient for each requried Legendre polynomial
void get_a_hat (legendre_info &basis, 
		particle_info &a)
{

for (int i = 0; i<basis.M; i++)
{
	basis.a_hat_n[i] = basis.A_n[i] / basis.N;
	basis.a_hat_m[i] = basis.A_m[i] / basis.N;
	basis.a_hat_n_m[i] = basis.A_n_m[i] / basis.N;

	basis.sigma_a_n_a_m[i] = basis.a_hat_n_m[i] - basis.a_hat_n[i] * 
	basis.a_hat_m[i];
}
}

//Solves for the orthogonality constant due to the phase space shift
float get_ortho_const(int n, 
		      legendre_info & basis)
{
	return (2.0*n+1.0)/(basis.max-basis.min);
}

// Solves for the current and the uncertainty due to each coefficient
void get_current (legendre_info &basis)
{
	for(int n=0; n<basis.M; n++)
	{
		basis.ortho_const_n[n] = get_ortho_const(n,basis);
		basis.current[n] = basis.a_hat_n[n] * basis.ortho_const_n[n];

		for(int m=0;m<basis.M;m++)
		{
			basis.ortho_const_m[m] = get_ortho_const(m,basis);
			basis.current_unc[m] += basis.ortho_const_n[n] *
			basis.ortho_const_m[m] * basis.sigma_a_n_a_m[m];
		}
	}
}
/*
void get_current (legendre_info &basis)
{

	for(int n=0; n<basis.M, n++)
	{
	basis.ortho_const[i] = (2*n+1)/(basis.max-basis.min);
	basis.current[i] = basis.a_hat_n[i] * basis.ortho_const[i] * Pn(n,
	}

}*/
//A new line

