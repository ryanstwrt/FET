//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   Shift/mc_tallies/FETs/legendre_solver.cpp
 * \author Ryan H. Stewart
 * \date   Wed June 27 05:08:30 2017
 * \brief  Function for solving legendre polynomial.
 */
//---------------------------------------------------------------------------//


#include<iostream>
#include "FET.hh"

using namespace std;

//for n=0;
float P0(float x)
{
	return 1.0;
}

//for n=1
float P1(float x)
{
	return x;
}

//for n
float Pn(int n, float x)
{
	if (n==0)
	{
	return P0(x);
	}

	else if (n==1)
	{
	return P1(x);
	}
//calls out the special case of x=1
	if (x == 1.0)
	{
	return 1.0;
	}

	if (x == -1.0)
	{
	return ((n % 2 ==0) ? 1.0 : -1.0);
	}
	
	if ((x==0) && (n % 2))
	{
	return 0.0;
	}

	float pn;
	pn = ((2*n-1)*x*Pn(n-1,x)-(n-1)*Pn(n-2,x))/n;

	return pn;	
}
	

