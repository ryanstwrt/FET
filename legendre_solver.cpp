//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   Shift/mc_tallies/FETs/legendre_solver.cpp
 * \author Ryan H. Stewart
 * \date   Wed June 27 05:08:30 2017
 * \brief  Function for solving legendre polynomial.
 */
//---------------------------------------------------------------------------//


#include "FET.hh"

using namespace std;

//for the special case of n=0;
double P0(double x)
{
	return 1;
}

//for the special case of n=1

double P1(double x)
{
	return x;
}

//for any Legendre polynomial in the nth order
double Pn(int n, 
	  double x)
{
	if (n==0)
	{
	return P0(x);
	}
	else if (n==1)
	{
	return P1(x);
	}
	else
	{
	double pn;
	pn = ((2*n-1)*x*Pn(n-1,x)-(n-1)*Pn(n-2,x))/n;
	return pn;
	}
}
