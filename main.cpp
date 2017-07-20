

#ifndef FET_H
#define FET_H
#include"FET.hh"
#endif

int main ()
{

tally_info tally;
legendre_info basis;
particle_info a;

tally.surface_index = 1;

//test test test
initialize_tally_info (tally, basis);
std::cout<<"test"<<std::endl;
surface_eval (basis, a, tally);


for(int n=0; n<basis.M; n++)
{	

std::cout<<basis.current[n]<<"*P_"<<n<<"(x) +/- "<<basis.var_a_n[n]<<std::endl;

}

std::cout<<std::endl;

}
