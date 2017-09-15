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
void FET_solver::get_current (legendre_info &basis,
	    	tally_info &tally, 
		std::size_t poly_terms,
		std::size_t N)
{
    std::vector<double> ortho_const(poly_terms, 0.0);
    std::vector<double> var_a_n_x(poly_terms, 0.0);
    std::vector<double> var_a_n_y(poly_terms, 0.0);

    std::vector<double> var_b_n_x(poly_terms, 0.0);
    std::vector<double> var_b_n_y(poly_terms, 0.0);
    std::vector<double> var_b_n_z(poly_terms, 0.0);

    std::vector<double> var_b(poly_terms, 0.0);

    for (int m = 0; m<poly_terms; m++)
    {
//Solve for the variance and coefficient vectors for surface tallies
        var_a_n_x[m] = (basis.A_n_x[m+poly_terms] - (1.0/N)*(basis.A_n_x[m] * basis.A_n_x[m])) * 1.0/(N*(N-1.0));
        var_a_n_y[m] = (basis.A_n_y[m+poly_terms] - (1.0/N)*(basis.A_n_y[m] * basis.A_n_y[m])) * 1.0/(N*(N-1.0));
        basis.A_n_x[m] *= (basis.x_basis[basis.surface_index+1]-basis.x_basis[basis.surface_index]) / N;
        basis.A_n_y[m] *= (basis.y_basis[basis.surface_index+1]-basis.y_basis[basis.surface_index]) / N;

        ortho_const[m] = (2.0*m+1.0)/2.0;
        tally.current_matrix_x[m] = basis.A_n_x[m] * ortho_const[m];
	tally.current_matrix_y[m] = basis.A_n_y[m] * ortho_const[m];
        tally.current_unc_matrix[m] = std::sqrt(var_a_n_x[m]);
        tally.current_R_sqr_value[m] = (var_a_n_x[m] * ortho_const[m]) / std::pow(basis.A_n_x[m],2.0);
	tally.current_total_unc += std::sqrt(pow(var_a_n_x[m],2) + pow(var_a_n_y[m],2));

//Solve for the variacne and the coefficient vectors for collision (volume) tallies


        var_b_n_x[m] = (basis.B_n_x[m+poly_terms] - (1.0/N)*(basis.B_n_x[m] * basis.B_n_x[m])) * 1.0/(basis.n_counter[0]*(basis.n_counter[0]-1.0));
        var_b_n_y[m] = (basis.B_n_y[m+poly_terms] - (1.0/N)*(basis.B_n_y[m] * basis.B_n_y[m])) * 1.0/(basis.n_counter[0]*(basis.n_counter[0]-1.0));
        var_b_n_z[m] = (basis.B_n_y[m+poly_terms] - (1.0/N)*(basis.B_n_z[m] * basis.B_n_z[m])) * 1.0/(basis.n_counter[0]*(basis.n_counter[0]-1.0));


//	var_b[m] = (basis.B[m+poly_terms] - (1.0/N)*(basis.B[m] * basis.B_n_x[m])) * 1.0/(basis.n_counter[0]*(basis.n_counter[0]-1.0));
   for(int n=0; n<poly_terms; ++n)
    {
      for(int i=0; i<poly_terms; ++i)
      {
	basis.B[m][n][i] *= (basis.x_basis[basis.surface_index+1]-basis.x_basis[basis.surface_index]) / basis.n_counter[0];
	tally.flux_matrix[m][n][i] = basis.B[m][n][i] * ortho_const[m] * ortho_const[n] * ortho_const[i];
      }
    }
//	tally.flux_unc[m] = std::sqrt(fabs(var_b[m]));
//	tally.flux_R_sqr[m] = (var_b[m] * ortho_const[m]) / std::pow(basis.B[m],2.0);
    }
}

