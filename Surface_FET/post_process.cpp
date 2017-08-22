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
    std::vector<double> ortho_const(poly_terms, 0.0);
    std::vector<double> var_a_n_x(poly_terms, 0.0);
    std::vector<double> var_a_n_y(poly_terms, 0.0);


    for (int m = 0; m<poly_terms; m++)
    {
        var_a_n_x[m] = (basis.A_m_x[m] - (1.0/N)*(basis.A_n_x[m] * basis.A_n_x[m])) * 1.0/(N*(N-1.0));
        var_a_n_y[m] = (basis.A_m_y[m] - (1.0/N)*(basis.A_n_y[m] * basis.A_n_y[m])) * 1.0/(N*(N-1.0));
        basis.A_n_x[m] *= (basis.x_basis[basis.surface_index+1]-basis.x_basis[basis.surface_index]) / N;
        basis.A_n_y[m] *= (basis.y_basis[basis.surface_index+1]-basis.y_basis[basis.surface_index]) / N;

        ortho_const[m] = (2.0*m+1.0)/2.0;
        tally.current_matrix_x[m] = basis.A_n_x[m] * ortho_const[m];
	tally.current_matrix_y[m] = basis.A_n_y[m] * ortho_const[m];
        tally.unc_matrix[m] = std::sqrt(var_a_n_x[m]);
        tally.R_sqr_value[m] = (var_a_n_x[m] * ortho_const[m]) / std::pow(basis.A_n_x[m],2.0);
    }
}

