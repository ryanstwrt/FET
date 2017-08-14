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
void get_current (legendre_info &basis,
	    	tally_info &tally, 
		std::size_t poly_terms,
		std::size_t N)
{
    for (int m = 0; m<poly_terms; m++)
    {
        basis.var_a_n[m] = (basis.A_m[m] - (1.0/N)*(basis.A_n[m] * basis.A_n[m])) * 1.0/(N*(N-1.0));
        basis.A_n[m] *= (basis.max-basis.min) / N;

        basis.ortho_const[m] = (2.0*m+1.0)/2.0;
        tally.current_matrix[m] = basis.A_n[m] * basis.ortho_const[m];
        tally.unc_matrix[m] = std::sqrt(basis.var_a_n[m]);
        tally.R_sqr_value[m] = (basis.var_a_n[m] * basis.ortho_const[m]) / std::pow(basis.A_n[m],2.0);
    }
}

