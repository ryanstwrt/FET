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
	}

	for (int n=0; n < basis.N; n++)
	{
		for(int m=0; m<basis.M; m++)
		{
		basis.sigma_a_n_a_m[m][n] = 
		basis.A_n_m[m][n]/basis.N - (basis.a_hat_n[n] *
		basis.a_hat_m[m]);
		}
	}
}

//Solves for the orthogonality constant due to the phase space shift
float legendre_info::get_ortho_const(int n, 
		      legendre_info & basis)
{
	return (2.0*n+1.0)/(basis.max-basis.min);
}

// Solves for the current and the uncertainty due to each coefficient

void legendre_info::get_current (legendre_info &basis)
{
	for(int m=0; m<basis.M; m++)
	{
		basis.ortho_const_n[m] = basis.get_ortho_const(m,basis);
		basis.current[m] = basis.a_hat_n[m];
	}
	
	for(int m=0;m<basis.M;m++)
	{
		for(int n=0;n<basis.N;n++)
		{
			basis.current_unc[m] += basis.sigma_a_n_a_m[m][n];
		}
		basis.current_unc[m] *= basis.N/(basis.N-1);
	}
}
