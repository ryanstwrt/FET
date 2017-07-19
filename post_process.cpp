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
	std::vector<double> term1;
	for (int i = 0; i<basis.M; i++)
	{
		basis.a_hat_n[i] = basis.A_n[i]/basis.N;
		basis.a_hat_m[i] = basis.A_m[i]/basis.N;
	}

	for (int n=0; n < basis.N; n++)
	{
		for(int m=0; m<basis.M; m++)
		{
			if(n==0)
			{
			term1.push_back(0);
			}
			term1[m] += std::pow(basis.A_n_m[m][n],2);
			std::cout<<basis.A_n_m[m][n]<<"  ";
/*			basis.sigma_a_n_a_m[m][n] = basis.A_n_m[m][n]/basis.N - (basis.a_hat_n[n] * basis.a_hat_m[m]);
			if(n==basis.N-1)
			{
				basis.var_a_n[m] = pow(basis.var_a_n[m],2) / pow(basis.A_n[m],2);
				basis.var_a_n[m] = std::sqrt(std::abs(basis.var_a_n[m]));
			}
*/		}
std::cout<<std::endl;
	}
		for(int m=0; m<basis.M; m++)
		{
			basis.var_a_n[m] = (term1[m]  - std::pow(basis.a_hat_n[m],2)/basis.N)/(basis.N*(basis.N-1));
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
	for(int m=0; m<basis.M; m++)
	{
		basis.ortho_const_n[m] = get_ortho_const(m,basis);
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
