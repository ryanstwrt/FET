//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   Shift/mc_tallies/FETs/surface_FET.cpp
 * \author Ryan H. Stewart
 * \date   Wed June 28 04:30:47 2017
 * \brief  Funciton for scaling to legendre space.
 */
//---------------------------------------------------------------------------//

#include"FET.hh"

//Solves for the Legendre coefficient and the associated uncertainty for each coefficient(currently broken)
void get_a_hat (legendre_info &basis, 
		particle_info &a,
	    	tally_info &tally)
{
	std::vector<double> term1;
	std::vector<double> term1a;
	for (int m = 0; m<basis.M; m++)
	{
		basis.a_hat_n[m] = basis.A_n[m]/basis.N;
		basis.a_hat_m[m] = basis.A_m[m]/basis.N;
		basis.a_hat_n1[m] = tally.surface_tallies[m][basis.n_counter][tally.surface_index]/basis.N;
		std::cout<<basis.a_hat_n[m]<<"  "<<basis.a_hat_n1[m]<<std::endl;
//		basis.a_hat_n1[m] = tally.surface_tallies[m][basis.n_counter][tally.surface_index]/basis.N;

	}

	for (int n=0; n < basis.N; n++)
	{
		for(int m=0; m<basis.M; m++)
		{
			if(n==0)
			{
			term1.push_back(0);
			term1a.push_back(0);
			}
			term1[m] += std::pow(basis.A_n_m[m][n],2);
			term1a[m] += std::pow(tally.surface_tallies[m][n][0],2);
//			term1a[m] += std::pow(tally.surface_tallies[m][n][tally.surface_index],2);
//			std::cout<<term1a[m]<<"  ";
//			std::cout<<term1a[m]<<std::endl;
/*			basis.sigma_a_n_a_m[m][n] = basis.A_n_m[m][n]/basis.N - (basis.a_hat_n[n] * basis.a_hat_m[m]);
			if(n==basis.N-1)
			{
				basis.var_a_n[m] = pow(basis.var_a_n[m],2) / pow(basis.A_n[m],2);
				basis.var_a_n[m] = std::sqrt(std::abs(basis.var_a_n[m]));
			}
*/		}
	}
		for(int m=0; m<basis.M; m++)
		{
			
			basis.var_a_n[m] = (term1[m]  - std::pow(basis.a_hat_n[m],2)/basis.N)/(basis.N*(basis.N-1));
			basis.var_a_n1[m] = (term1[m] - std::pow(basis.a_hat_n1[m],2)/basis.N)/(basis.N*(basis.N-1));
		}
}

//Solves for the orthogonality constant due to the phase space shift
float get_ortho_const(int n, 
		      legendre_info & basis)
{
	return (2.0*n+1.0)/(basis.max-basis.min);
}

// Solves for the current and the overall uncertainty

void get_current (legendre_info &basis,
	  	  tally_info &tally)
{
	
	for(int m=0; m<basis.M; m++)
	{
		basis.ortho_const_n[m] = get_ortho_const(m,basis);
		basis.current[m] = basis.a_hat_n[m];
		basis.current1[m] = basis.a_hat_n1[m];
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
