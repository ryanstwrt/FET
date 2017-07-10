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
#include"FET.hh"

void get_a_hat (legendre_info &basis, particle_info &a)
{

for (int i = 0; i<basis.M; i++)
{
basis.a_hat_n[i] = basis.A_n[i] / basis.N;
basis.a_hat_m[i] = basis.A_m[i] / basis.N;
basis.a_hat_n_m[i] = basis.A_n_m[i] / basis.N;

basis.sigma_a_n_a_m[i] = basis.a_hat_n_m[i] - basis.a_hat_n[i] * basis.a_hat_m[i];
}

}

void get_current (legendre_info &basis)
{

	for(int n=0; n<basis.M, n++)
	{
	basis.ortho_const[i] = (2*n)/(basis.max-basis.min)
	basis.current[i] = basis.a_hat_n[i] * basis.ortho_const[i] * Pn(n,
	}

}

