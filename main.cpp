

#ifndef FET_H
#define FET_H
#include"FET.hh"
#endif

int main ()
{

legendre_info basis;
particle_info a;

initalize (basis, a);
a.get_particle(a);


basis_eval(basis, a);


get_a_hat(basis, a);
basis.get_current(basis);

std::cout<<"The scaled current is described as: "<<std::endl;
for(int n=0; n<basis.M; n++)
{
	if(n==0)
	{
		std::cout<<basis.current[n]<<" +/- "<<basis.var_a_n[n]<<std::endl;
	}
	else
	{
		std::cout<<basis.current[n]<<"*x^"<<n<<" +/- "<<basis.var_a_n[n]<<std::endl;
	}
}

std::cout<<std::endl;
}
