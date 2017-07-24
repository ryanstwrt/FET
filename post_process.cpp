//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   Shift/mc_tallies/FETs/surface_FET.cpp
 * \author Ryan H. Stewart
 * \date   Wed June 28 04:30:47 2017
 * \brief  Funciton for scaling to legendre space.
 */
//---------------------------------------------------------------------------//

#include"FET.hh"

//Solves for the Legendre coefficient and the associated uncertainty for each coefficient
void get_a_hat (legendre_info &basis, 
		particle_info &a,
	    	tally_info &tally)
{

	for (int m = 0; m<basis.M; m++)
	{
		basis.a_hat_n[m] = basis.A_n[m]/basis.N;
		basis.a_hat_m[m] = basis.A_m[m]/basis.N;
		basis.a_hat_n1[m] = tally.surface_tallies[m][basis.n_counter][tally.surface_index]/basis.N;
	}

}

void get_uncertainty (legendre_info &basis, 
		      particle_info &a,
		      tally_info &tally)
{

	std::vector<double> term1;
	std::vector<double> term1a;
	term1a.resize(basis.M);

//Get the first term in the uncertainty	
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
			term1a[m] += std::pow(tally.surface_tallies[m][n][tally.surface_index],2);
		}
	}

//Get the uncertainty for each coefficient
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
	          particle_info &a,
	  	  tally_info &tally)
{

	for(int s=0; s<tally.num_surfaces; s++)
	{
		tally.surface_index = s;	
		get_a_hat(basis, a, tally);
		get_uncertainty(basis, a, tally);
	
		for(int m=0; m<basis.M; m++)
		{
			basis.ortho_const_n[m] = get_ortho_const(m,basis);
			basis.total_current[m] = basis.a_hat_n[m];
			basis.current1[m] = basis.a_hat_n1[m];
			tally.current_matrix[m][s] = basis.a_hat_n1[m];
			tally.unc_matrix[m][s] = basis.var_a_n1[m];
			//std::cout<<basis.total_current[m]<<"  "<<basis.current1[m]<<"  "<<tally.current_matrix[m][s]<<std::endl;
			std::cout<<tally.current_matrix[m][s]<<" +/- "<<tally.unc_matrix[m][s]<<std::endl;
		}
		std::cout<<std::endl;
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
