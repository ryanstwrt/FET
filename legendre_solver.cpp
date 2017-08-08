//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   Shift/mc_tallies/FETs/legendre_solver.cpp
 * \author Ryan H. Stewart
 * \date   Wed June 27 05:08:30 2017
 * \brief  Function for solving legendre polynomial.
 */
//---------------------------------------------------------------------------//


#include "FET.hh"



//for the special case of n=0;
//verified 7/18/17
double P0(double x)
{
    return 1;
}

//for the special case of n=1
//verified 7/18/17
double P1(double x)
{
    return x;
}

double P2(double x)
{
    
}

double P3(double x)
{

}

double P4(double x)
{

}

double P5(double x)
{

}

double P6(double x)
{

}

double P7(double x)
{

}

double P8(double x)
{

}

double P9(double x)
{

}

double P10(double x)
{

}

double P11(double x)
{

}

double P12(double x)
{

}
//for any Legendre polynomial in the nth order
//verified 7/18/17
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
