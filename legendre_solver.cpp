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
const double x2 = x * x;
double pn;
/*    if (n<=12)
    {
    switch(n)
    {
	case 12:
	pn = ((((((676039 * x2 - 1939938) * x2 + 2078505) * x2 - 1021020) * x2 + 225225) * x2 - 18018) * x2 + 231) / 1024;
	case 11:
	pn = (((((88179 * x2 - 230945) * x2 + 218790) * x2 - 90090) * x2 + 15015) * x2 - 693) * x / 256;
	case 10:
	pn = (((((46189 * x2 - 109395) * x2 + 90090) * x2 - 30030) * x2 + 3465) * x2 - 63) / 256;
	case 9:
	pn = ((((12155 * x2 - 25740) * x2 + 18018) * x2 - 4620) * x2 + 315) * x / 128;
	case 8:
	pn = ((((6435 * x2 - 12012) * x2 + 6930) * x2 - 1260) * x2 + 35) / 128;
	case 7:
	pn = (((429 * x2 - 693) * x2 + 315) * x2 - 35) * x / 16;
	case 6:
	pn = (((231 * x2 - 315) * x2 + 105) * x2 - 5) / 16;
	case 5:
	pn = ((35 * x2 - 30) * x2 + 3) / 8; 
	case 4:
	pn = ((35 * x2 - 30) * x2 + 3) / 8;
	case 3:
	pn = (5 * x2 - 3) * x / 2;
	case 2:
	pn = (3 * x2 - 1) / 2;
	case 1:
	pn = x;
	case 0:
	pn = 1.0;
    }
    else*/
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
