//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   Shift/mc_tallies/FETs/surface_FET.cpp
 * \author Ryan H. Stewart
 * \date   Wed June 28 04:30:47 2017
 * \brief  Funciton for scaling to legendre space.
 */
//---------------------------------------------------------------------------//

#include"FET.hh"

//Solves for each Legendre coefficient
void get_legendre_coefficient (legendre_info &basis, 
		particle_info &a,
	    	tally_info &tally)
{
    for(int n=0; n<basis.N; n++)
    {
        for (int m = 0; m<basis.M; m++)
	{
	    if(n==0)
	    {
		basis.A_n[m] = 0;
	        basis.unc_term1[m] = 0;
	    }

	    basis.A_n[m] += std::fabs(tally.surface_tallies[m][n][tally.surface_index]);
	    basis.unc_term1[m] += std::pow(tally.surface_tallies[m][n][tally.surface_index],2);
//	    std::cout<<tally.surface_tallies[m][n][tally.surface_index]<<"  "<<basis.A_n[m]<<"  "<<basis.unc_term1[m]<<std::endl;
	}
	//std::cout<<std::endl;
    }

    for (int m = 0; m<basis.M; m++)
    {
        basis.a_hat_n[m] = basis.A_n[m]/basis.N;
        basis.var_a_n[m] = (basis.unc_term1[m] - (std::pow(basis.A_n[m],2)/basis.N) )/ (basis.N*(basis.N-1));
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
    for(tally.surface_index=0; tally.surface_index<tally.num_surfaces;
    tally.surface_index++)
    {
	get_legendre_coefficient(basis, a, tally);
	
	for(int m=0; m<basis.M; m++)
	{
	    basis.ortho_const[m] = get_ortho_const(m,basis);
	    tally.coefficient_matrix[m][tally.surface_index] = basis.a_hat_n[m];
	    tally.current_matrix[m][tally.surface_index] = basis.a_hat_n[m] * basis.ortho_const[m];
	    tally.unc_matrix[m][tally.surface_index] = std::sqrt(fabs(basis.var_a_n[m]));//  / basis.a_hat_n[m];
 	    tally.R_sqr_value[m][tally.surface_index] = (std::fabs(basis.var_a_n[m]) * basis.ortho_const[m]) / std::pow(basis.a_hat_n[m],2);

	    //std::cout<<basis.ortho_const[m]<<std::endl;
	    std::cout<<tally.current_matrix[m][tally.surface_index]<<"  ";
	    std::cout<<basis.unc_term1[m]<<"  ";
   	    std::cout<<tally.coefficient_matrix[m][tally.surface_index]<<std::endl;
std::cout<<std::endl;
	}
std::cout<<std::endl;
    }
}
