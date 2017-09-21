/* 
 * File:   Distribution.h
 * Author: brycen
 *
 * Created on August 2, 2017, 12:34 PM
 */

#ifndef DISTRIBUTION_H
#define DISTRIBUTION_H

class Distribution
{
public:
  Distribution() = default;
  virtual ~Distribution() = default;
  
  double operator() (double x) const;
};

#endif /* DISTRIBUTION_H */

