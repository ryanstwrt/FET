//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   Shift/mc_tallies/FET.cc
 * \author Ryan H. Stewart
 * \date   Wed June 28 04:30:47 2017
 * \mod    Thur Sept 25 09:48:00 2017
 * \brief  Collection of FET functions for solving
 */
//---------------------------------------------------------------------------//

#include"FET.hh"

using namespace std;

void FET_solver::initializer (initial_info &info)
{

info.poly_order = 9;
info.poly_terms = info.poly_order + 1;
info.num_tallies = 1;


}
//---------------------------------------------------------------------------//
/*!
 * Below solves for the contribution of an individual particle to the tally.
 *
 * 
 */
//---------------------------------------------------------------------------//
//Scales the reactor phase space down to Legendre phase space [-1,1]
//Verified 7/18/17
double FET_solver::scale (double position, 
	     std::vector<double> dimension_basis)
{
    if(position < dimension_basis[0] || position > dimension_basis[1])
    {
	return 2;
    }
    else
    {
	return (2 * ((position - dimension_basis[0])/(dimension_basis[1] - dimension_basis[0])) - 1);

    }
}

//Solves the legendre polynomial for a particles nth contribution to an estimate of the coefficient i.e. this will happen multiple times for each particle until the particle dies
void FET_solver::surface_eval (legendre_info &basis, 
		   particle_info &a, std::size_t poly_terms)
{
    double x_tild = scale(a.x, basis.x_basis);
    double y_tild = scale(a.y, basis.y_basis);

  if(x_tild != 2 || y_tild != 2)
  {
    std::vector<double> P_n_x = basis.Pn(poly_terms, x_tild);
    std::vector<double> P_n_y = basis.Pn(poly_terms, y_tild);

//Sums the contribution to the current of each collision from one particle

    for(int m=0; m<poly_terms; ++m)
    {
	for(int n=0; n<poly_terms; ++n)
	{
          basis.a[m][n] += a.wt * P_n_x[m] * P_n_y[n];
        }
    }
  }
}

//Particle death triggers eval and generates an estimate for the coefficient for the nth particle.
void FET_solver::surface_eval2 (legendre_info &basis, 
		   particle_info &a, std::size_t poly_terms)
{
//Sums the contribution to the current of each particle
    for(int m=0; m<poly_terms; ++m)
    {
	for(int n=0; n<poly_terms; ++n)
	{
	  basis.A[m][n] += basis.a[m][n];
	  basis.A_unc[m][n] += pow(basis.a[m][n],2);
	  basis.a[m][n] = 0;
	}
    }

   basis.n_counter[0]++;
}

//Solves the legendre polynomial for a particles nth contribution to an estimate of the coefficient i.e. this will happen multiple times for each particle until the particle dies
//TO DO: Combine collision eval with an if statement (if p.alive == 1)
void FET_solver::collision_eval (legendre_info &basis, 
		   particle_info &a, std::size_t poly_terms)
{
    double x_tild = scale(a.x, basis.x_basis);
    double y_tild = scale(a.y, basis.y_basis);
    double z_tild = scale(a.z, basis.z_basis);
  
  if(x_tild != 2 || y_tild != 2 || z_tild != 2 )
  {
    std::vector<double> P_n_x = basis.Pn(poly_terms, x_tild);
    std::vector<double> P_n_y = basis.Pn(poly_terms, y_tild);
    std::vector<double> P_n_z = basis.Pn(poly_terms, z_tild);

//Sums the contribution to the flux of each collision from one particle
    for(int m=0; m<poly_terms; ++m)
    {
      for(int n=0; n<poly_terms; ++n)
      {
        for(int i=0; i<poly_terms; ++i)
        {
          basis.b[m][n][i] += a.wt * P_n_x[m] * P_n_y[n] * P_n_z[i] / a.xs_tot;
        }
      }
    }
  }
}

//Particle death triggers eval and generates an estimate for the coefficient for the nth particle.
void FET_solver::collision_eval2(legendre_info &basis, 
		   particle_info &a, std::size_t poly_terms)
{
//Sums the contribution of the flux of each particle
    for(int m=0; m<poly_terms; ++m)
    {
      for(int n=0; n<poly_terms; ++n)
      {
        for(int i=0; i<poly_terms; ++i)
        {
	  basis.B[m][n][i] += basis.b[m][n][i];
	  basis.B_unc[m][n][i] += pow(basis.b[m][n][i],2);
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
//---------------------------------------------------------------------------//


void FET_solver::get_current (legendre_info &basis,
	    	tally_info &tally, 
		std::size_t poly_terms)
{
//Dummy variables for solving for the current
    std::vector<std::vector<double> >   var_a;
    std::vector<std::vector<std::vector<double> > >  var_b;
    std::vector<std::vector<std::vector<double> > >   ortho_const;

//Initialize the dummy variables to the same dimensions as the current/flux
   var_b.resize(poly_terms);
   var_a.resize(poly_terms);
   ortho_const.resize(poly_terms);
   for(int m=0; m<poly_terms; ++m)
   {
     var_b[m].resize(poly_terms);
     var_a[m].resize(poly_terms);
     ortho_const[m].resize(poly_terms);
     for(int n=0; n<poly_terms; ++n)
     {
       ortho_const[m][n].resize(poly_terms);
       var_b[m][n].resize(poly_terms);
       var_a[m][n] = 0;
       for(int i=0; i<poly_terms; ++i)
       {
         var_b[m][n][i] = 0;
	 ortho_const[m][n][i] = (2*m-1) * (2*n-1)* (2*i-1) / ((basis.x_basis[1]-basis.x_basis[0])*(basis.y_basis[1]-basis.y_basis[0])*(basis.z_basis[1]-basis.z_basis[0]));
	 
       }
     }
   }

//Solves for the final coefficient, followed by the current, the uncertainty, and finally the R^2 value
  for (int m=0; m<poly_terms; m++)
  {
    for(int n=0; n<poly_terms; ++n)
    {
      basis.A[m][n] *= ((basis.x_basis[0]-basis.x_basis[1])*(basis.y_basis[0]-basis.y_basis[1])) / basis.n_counter[0];
      tally.current_matrix[m][n] = basis.A[m][n] * ortho_const[m][n][0];

      var_a[m][n] = (basis.A_unc[m][n] - (1.0/basis.n_counter[0]) * std::pow(basis.A[m][n],2) ) * 1.0 / (basis.n_counter[0]*(basis.n_counter[0]-1.0));
      tally.current_unc_matrix[m][n] = std::sqrt(fabs(var_a[m][n]));
      tally.current_R_matrix[m][n] = (var_a[m][n] * ortho_const[m][n][0] ) / std::pow(basis.A[m][n],2.0);
    }	
  }

//Solves for the final coefficient, followed by the flux, the uncertainty, and finally the R^2 value
  for(int m=0; m<poly_terms; ++m)
  {
    for(int n=0; n<poly_terms; ++n)
    {
      for(int i=0; i<poly_terms; ++i)
      {
	basis.B[m][n][i] *= 1 / basis.n_counter[0];
        var_b[m][n][i] = (basis.B_unc[m][n][i]-(1.0/basis.n_counter[0])*std::pow(basis.B[m][n][i],2))*1.0/(basis.n_counter[0]*(basis.n_counter[0]-1.0));


	tally.flux_R_matrix[m][n][i] = (var_b[m][n][i])  * ortho_const[m][n][i]/ std::pow(basis.B[m][n][i],2.0);
	tally.flux_matrix[m][n][i] = basis.B[m][n][i] * ortho_const[m][n][i]; 
	tally.flux_unc_matrix[m][n][i] = std::sqrt(fabs(var_b[m][n][i]));


      }
    }
  }
}

//Removes any coefficient with an R^2 value greater than 10. And lets the user know how many coefficients had R^2 values greater than 1 and greater than 10
void FET_solver::cleanup (tally_info &tally, 
		std::size_t poly_terms)
{

    tally.R_great_10 = 0;
    tally.R_great_1 = 0;
    tally.total_coeff = 0;

for (int m=0; m < poly_terms; ++m)
{
   for(int n=0; n<poly_terms; ++n)
    {
      for(int i=0; i<poly_terms; ++i)
      {
	if(tally.flux_R_matrix[m][n][i] >= 1)
	{
	    tally.R_great_1++;
	 if(tally.flux_R_matrix[m][n][i] >= 10)
	  {
	    tally.R_great_10++;
	    std::cout<<"P("<<m<<")("<<n<<")("<<i<<") = "<<tally.flux_matrix[m][n][i]<<" +/- "<<tally.flux_unc_matrix[m][n][i]<<" w/ "<<tally.flux_R_matrix[m][n][i]<<std::endl;
	    tally.flux_matrix[m][n][i] = 0;

	  }
	}
	  tally.total_coeff++;
      }
    }
    std::cout<<std::endl;
  }
}

//---------------------------------------------------------------------------//
/*!
 * Below solves the Legendre Polynomials for the highest order term and saves all of the lower order terms in a vector.
 *
 * 
 */
//---------------------------------------------------------------------------//

std::vector<double> legendre_info::Pn(std::size_t poly_terms, 
		       double x)
{
const double x2 = x * x;

//Switch case to solve for the first twelve cases of the Legendre polynomials
//This decreases the time requried for calculations
    switch(poly_terms)
    {
	default:
	case 12:
	save(12, ((((((676039 * x2 - 1939938) * x2 + 2078505) * x2 - 1021020) * x2 + 225225) * x2 - 18018) * x2 + 231) / 1024);
	case 11:
	save(11, (((((88179 * x2 - 230945) * x2 + 218790) * x2 - 90090) * x2 + 15015) * x2 - 693) * x / 256);
	case 10:
	save(10, (((((46189 * x2 - 109395) * x2 + 90090) * x2 - 30030) * x2 + 3465) * x2 - 63) / 256);
	case 9:
	save(9, ((((12155 * x2 - 25740) * x2 + 18018) * x2 - 4620) * x2 + 315) * x / 128);
	case 8:
	save(8, ((((6435 * x2 - 12012) * x2 + 6930) * x2 - 1260) * x2 + 35) / 128);
	case 7:
	save(7, (((429 * x2 - 693) * x2 + 315) * x2 - 35) * x / 16);
	case 6:
	save(6, (((231 * x2 - 315) * x2 + 105) * x2 - 5) / 16);
	case 5:
	save(5, ((63 * x2 - 70) * x2 + 15)*x / 8); 
	case 4:
	save(4, ((35 * x2 - 30) * x2 + 3) / 8);
	case 3:
	save(3, (5 * x2 - 3) * x / 2);
	case 2:
	save(2, (3 * x2 - 1) / 2);
	case 1:
	save(1, x);
	case 0:
	save(0, 1.0);
    }
//Solves for any case above 12 and saves it in the polynomial vector
	for (int n = 13; n < poly_terms; ++n)
	{
	save(n, ( ( 2 * n - 1 ) * x * load( n - 1 ) - ( n - 1 ) * load( n - 2 ) ) / n );
	}

    return P_n;
}

//Serves as a loader and saver for the polynomial vector
double
legendre_info::load(std::size_t index) const
{
  return P_n[index];
}

void
legendre_info::save(std::size_t index, double value)
{
  P_n[index] = value;
}

