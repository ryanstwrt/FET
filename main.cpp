
#include<iostream>
#include<vector>
#include<cstdlib>
#include<cmath>
#include <cstdlib>
#include <ctime>
#include <iomanip>
#include <random>
#include <time.h>

#include"FET.hh"
#include "Distribution.hh"

int main (int argc, char** argv)
{
clock_t tStart = clock ();

const std::size_t poly_order = 10;
const std::size_t poly_terms = poly_order + 1;
const std::size_t N = 1e7;

Distribution fluxshape;

tally_info tally;
legendre_info basis;
particle_info a;
tally.surface_index=0;

initialize_tally_info (tally, basis);

for(int n = 0; n<basis.N; n++)
{
	a.get_particle(basis, a);

	surface_eval (basis, a, tally);
}

get_current(basis, a, tally);

for(int s=0; s<tally.num_surfaces; s++)
{	
	std::cout<<"Surface Number "<<s<<std::endl;
	for(int m=0; m<basis.M; m++)
	{
		std::cout<<tally.coefficient_matrix[m][s]<<" * P_" <<m<<"(x) +/- "<<tally.unc_matrix[m][s]<<"    ";
		std::cout<<"With an R^2 value of "<<tally.R_sqr_value[m][s]<<std::endl;
	}

std::cout<<std::endl;
}

for(int cp = 0; cp <= 10; cp++)
{
 float x = double(cp) / 10 * 2 -1;
 float y = fluxshape(x);
 float sum = 0;

  for(int m = 0; m < basis.M; m++)
	sum += tally.coefficient_matrix[m][0] * Pn(m,x);
    std::cout.precision(2);
    std::cout << std::fixed;
    std::cout << "Position: " << x
              << "\t\tFlux Value: " << y;
    
    std::cout.precision(5);
    std::cout << "\t\tFE Value: " << sum
              << "\t\tDifference: " << std::scientific << y - sum << std::endl;
}
std::cout<<"Time taken: "<<(double)(clock() - tStart)/CLOCKS_PER_SEC<<std::endl;
}
