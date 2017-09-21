//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   Shift/mc_tallies/FET.cc
 * \author Ryan H. Stewart
 * \date   Wed June 28 04:30:47 2017
 * \mod    Thur Sept 20 04:53:00 2017
 * \brief  Collection of FET functions for solving
 */
//---------------------------------------------------------------------------//

#include"FET.hh"

using namespace std;

//Scales the original phase space down to Legendre phase space [-1,1]
//Verified 7/18/17
double FET_solver::scale (double x, 
	     legendre_info basis)
{
    if(x < basis.x_basis[basis.surface_index] || x > basis.x_basis[basis.surface_index+1])
    {
	std::cout<< "The domain does not encompass the entire range of the problem."<<endl;
    }
    else
    {
	return (2 * ((x - basis.x_basis[basis.surface_index])/(basis.x_basis[basis.surface_index+1] - basis.x_basis[basis.surface_index])) - 1);

    }
}

//Calculated the individual contribution from one particle, which can contribute multiple times depending on how many times it crosses the surface
//Verified 7/18/17
void FET_solver::surface_eval (legendre_info &basis, 
		   particle_info &a, std::size_t poly_terms)
{
    double temp;
    double x_tild = scale(a.b_weight, basis);
    double y_tild = scale(a.b_weight, basis);
    double ratio = a.k_particle;
    double ratio2 = a.k_particle;
    std::vector<double> a_n_x(poly_terms, 0.0);
    std::vector<double> a_n_y(poly_terms, 0.0);
    std::vector<double> P_n_x = basis.Pn(poly_terms, x_tild);
    std::vector<double> P_n_y = basis.Pn(poly_terms, y_tild);

    //calculate the legendre coefficients up to truncation value M for a_n and a_m for one particle
    //k is the number of times the particle crosses the specified surface, while m is the legendre coefficient
    for(int m=0; m<poly_terms; ++m)
    {
	for(int n=0; n<poly_terms; ++n)
	{
          basis.a[m][n] += a.b_weight * P_n_x[m] * P_n_y[n];
        }
    }

    for(int m=0; m<poly_terms; ++m)
    {
	for(int n=0; n<poly_terms; ++n)
	{
  	  temp =  basis.a[m][n];
	  basis.A[m][n] += temp;
	  basis.A_unc[m][n] += pow(temp,2);
	  basis.a[m][n] = 0;
	}
    }

   basis.n_counter[0]++;
}

void FET_solver::collision_eval (legendre_info &basis, 
		   particle_info &a, std::size_t poly_terms)
{

    double x_tild = scale(a.x, basis);
    double y_tild = scale(a.y, basis);
    double z_tild = scale(a.z, basis);
    std::vector<double> P_n_x = basis.Pn(poly_terms, x_tild);
    std::vector<double> P_n_y = basis.Pn(poly_terms, y_tild);
    std::vector<double> P_n_z = basis.Pn(poly_terms, z_tild);

//Calculates b with the coordinates from teh current neutron
  for(int m=0; m<poly_terms; ++m)
  {
    for(int n=0; n<poly_terms; ++n)
    {
      for(int i=0; i<poly_terms; ++i)
      {
        basis.b[m][n][i] += a.b_weight * P_n_x[m] * P_n_y[n] * P_n_z[i] / a.xs_tot;
      }
    }
  }

}

void FET_solver::collision_eval2(legendre_info &basis, 
		   particle_info &a, std::size_t poly_terms)
{
double temp;

    for(int m=0; m<poly_terms; ++m)
    {
      for(int n=0; n<poly_terms; ++n)
      {
        for(int i=0; i<poly_terms; ++i)
        {
  	  temp = basis.b[m][n][i];
	  basis.B[m][n][i] += temp;
	  basis.B_unc[m][n][i] += pow(temp,2);
	  basis.b[m][n][i] = 0;
        }
      }
   }

   basis.n_counter[0]++;
}

//---------------------------------------------------------------------------//
/*!
 * Below solves for the final current or the flux for the tally.
 *
 * 
 */

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

    std::vector<std::vector<double> >   var_a;
    std::vector<std::vector<std::vector<double> > >  var_b;

   var_b.resize(poly_terms);
   var_a.resize(poly_terms);
   for(int m=0; m<poly_terms; ++m)
   {
     ortho_const[m] = (2.0*m+1.0)/2.0;
     var_b[m].resize(poly_terms);
     var_a[m].resize(poly_terms);
     for(int n=0; n<poly_terms; ++n)
     {
       var_b[m][n].resize(poly_terms);
       var_a[m][n] = 0;
       for(int i=0; i<poly_terms; ++i)
       {
         var_b[m][n][i] = 0;

       }
     }
   }


  for (int m = 0; m<poly_terms; m++)
  {
    for(int n=0; n<poly_terms; ++n)
    {
      basis.A[m][n] *= (basis.x_basis[basis.surface_index+1]-basis.x_basis[basis.surface_index]) / basis.n_counter[0];
      tally.current_matrix[m][n] = basis.A[m][n] * ortho_const[m] * ortho_const[n];
      var_a[m][n] = (basis.A_unc[m][n] - (1.0/N) * std::pow(basis.A[m][n],2) ) * 1.0 / (basis.n_counter[0]*(basis.n_counter[0]-1.0));
      tally.current_unc_matrix[m][n] = std::sqrt(fabs(var_a[m][n]));
      tally.current_R_matrix[m][n] = (var_a[m][n] * ortho_const[m] * ortho_const[n] ) / std::pow(basis.A[m][n],2.0);
    }	
  }

  for(int m=0; m<poly_terms; ++m)
  {
    for(int n=0; n<poly_terms; ++n)
    {
      for(int i=0; i<poly_terms; ++i)
      {
	basis.B[m][n][i] *= (basis.x_basis[basis.surface_index+1]-basis.x_basis[basis.surface_index]) / basis.n_counter[0];
	tally.flux_matrix[m][n][i] = basis.B[m][n][i] * ortho_const[m] * ortho_const[n] * ortho_const[i];

        var_b[m][n][i] = (basis.B_unc[m][n][i] - (1.0/N) * std::pow(basis.B[m][n][i],2) ) * 1.0 / (basis.n_counter[0]*(basis.n_counter[0]-1.0)); 
	tally.flux_unc_matrix[m][n][i] = std::sqrt(fabs(var_b[m][n][i]));
	tally.flux_R_matrix[m][n][i] = (var_b[m][n][i] * ortho_const[m] * ortho_const[n] * ortho_const[i]) / std::pow(basis.B[m][n][i],2.0);
      }
    }
  }

}
