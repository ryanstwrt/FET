/* 
 * File:   Distribution.cpp
 * Author: brycen
 * 
 * Created on August 2, 2017, 12:34 PM
 */

#include "Distribution.hh"
#include "FET.hh"

double
Distribution::operator()(double x) const
{
  return 5 * x * x * x - x * x - 2 * x + 4;
}


