

#ifndef FET_H
#define FET_H
#include"FET.hh"
#endif

int main ()
{

legendre_info basis;
particle_info a;

basis.initalize (basis);
a.get_particle(basis, a);


basis_eval(basis, a);

get_a_hat(basis, a);
basis.get_current(basis);

std::cout<<"The scaled current is described as: "<<std::endl;
for(int n=0; n<basis.M; n++)
{
	std::cout<<basis.current[n]<<"*P_"<<n<<"(x) +/- "<<basis.var_a_n[n]<<std::endl;
}

std::cout<<std::endl;
}
