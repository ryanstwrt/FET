

#ifndef FET_H
#define FET_H
#include"FET.hh"
#endif

int main ()
{


tally_info tally;
legendre_info basis;
particle_info a;



//test test test

initialize_tally_info (tally, basis);


for(int n = 0; n<basis.N; n++)
{
a.get_particle(basis, a);
if(a.particle_surface <= 0.33)
{
	tally.surface_index=0;
}
else if (a.particle_surface >= 0.33 && a.particle_surface <=0.66)
{
	tally.surface_index=1;
}
else
{
	tally.surface_index=2;
}
	surface_eval (basis, a, tally);
	get_A (basis, a, tally);
//	std::cout<<tally.surface_index<<std::endl;

/*	for(int m=0; m<basis.M; m++)
	{
		std::cout<<a.a_n[m]<<"  ";
	}
	std::cout<<std::endl;*/
}

std::cout<<std::endl;

//get_a_hat(basis, a, tally);
get_current(basis, a, tally);
	

for(int n=0; n<basis.M; n++)
{	

std::cout<<basis.total_current[n]<<"*P_"<<n<<"(x) +/- "<<basis.var_a_n[n]<<std::endl;

}

std::cout<<std::endl;

}
