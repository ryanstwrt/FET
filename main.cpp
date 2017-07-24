

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
if(a.particle_surface <= 0.50)
{
	tally.surface_index=0;
}
else
{
	tally.surface_index=1;
}
	surface_eval (basis, a, tally);
	get_A (basis, a, tally);
	std::cout<<tally.surface_index<<std::endl;

	for(int m=0; m<basis.M; m++)
	{
		std::cout<<a.a_n[m]<<"  ";
	}
	std::cout<<std::endl;
}

std::cout<<std::endl;
tally.surface_index=0;
get_a_hat(basis, a, tally);
get_current(basis, tally);

for(int m=0; m<basis.M; m++)
{
	for(int s=0; s<tally.surface_index; s++)
	{
	std::cout<<tally.surface_tallies[m][basis.n_counter][s]<<"  ";
	}
std::cout<<std::endl;
}
	

for(int n=0; n<basis.M; n++)
{	

std::cout<<basis.current[n]<<"*P_"<<n<<"(x) +/- "<<basis.var_a_n[n]<<std::endl;

}

std::cout<<std::endl;

}
