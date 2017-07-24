

#ifndef FET_H
#define FET_H
#include"FET.hh"
#endif

int main (int argc, char *argv[])
{


tally_info tally;
legendre_info basis;
particle_info a;

tally.surface_index = 1;

//test test test
initialize_tally_info (tally, basis);


for(int n = 0; n<basis.N; n++)
{
	a.get_particle(basis, a);
	surface_eval (basis, a, tally);
	get_A (basis, a, tally);
}

get_a_hat(basis, a, tally);
get_current(basis, tally);

for(int n=0; n<basis.M; n++)
{	

std::cout<<basis.current[n]<<"*P_"<<n<<"(x) +/- "<<basis.var_a_n[n]<<std::endl;

}

std::cout<<std::endl;

}
