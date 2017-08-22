//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   Shift/mc_tallies/FETs/legendre_solver.cpp
 * \author Ryan H. Stewart
 * \date   Wed June 27 05:08:30 2017
 * \brief  Function for solving legendre polynomial.
 */
//---------------------------------------------------------------------------//


#include "FET.hh"

//for any Legendre polynomial in the nth order
//verified 7/18/17
std::vector<double> legendre_info::Pn(std::size_t poly_terms, 
		       double x)
{
const double x2 = x * x;

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
	save(5, ((35 * x2 - 30) * x2 + 3) / 8); 
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
	for (int n = 13; n < poly_terms; ++n)
	{
	save(n, ( ( 2 * n - 1 ) * x * load( n - 1 ) - ( n - 1 ) * load( n - 2 ) ) / n );
	}

    return P_n;
}

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
