

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
	else if (a.particle_surface <=0.66)
	{
		tally.surface_index=1;
	}
	else
	{
		tally.surface_index=2;
	}
	
	surface_eval (basis, a, tally);
}

get_current(basis, a, tally);

for(int s=0; s<tally.num_surfaces; s++)
{	
	std::cout<<"Surface Number "<<s<<std::endl;
	for(int m=0; m<basis.M; m++)
	{
		float test;
//		test = rescale(, basis);
		float test1;
//		test1 = rescale();
		std::cout<<tally.coefficient_matrix[m][s]<<" * P_" <<m<<"(x) +/- "<<tally.unc_matrix[m][s]<<"    ";
//		std::cout<<basis.total_current[m]<<"*P_"<<m<<"(x) +/- "<<basis.var_a_n[m]<<std::endl;
		std::cout<<"With an R^2 value of "<<tally.R_sqr_value[m][s]<<std::endl;
	}

std::cout<<std::endl;
}
}
