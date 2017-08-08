
#include<iostream>
#include<fstream>
#include<vector>
#include<cstdlib>
#include<cmath>
#include <cstdlib>
#include <ctime>
#include <iomanip>
#include <random>
#include <time.h>
using namespace  std;

#include"FET.hh"
#include "Distribution.hh"

int main (int argc, char** argv)
{
ofstream myfile;
myfile.open ("time.txt", ios::in | ios::app);
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
		std::cout<<tally.coefficient_matrix[m]<<" * P_" <<m<<"(x) +/- "<<tally.unc_matrix[m]<<"    ";
		std::cout<<"With an R^2 value of "<<tally.R_sqr_value[m]<<std::endl;
	}

std::cout<<std::endl;
}

for(int cp = 0; cp <= 10; cp++)
{
 float x = double(cp) / 10 * 2 -1;
 float y = fluxshape(x);
 float sum = 0;

  for(int m = 0; m < basis.M; m++)
	sum += tally.coefficient_matrix[m] * Pn(m,x);
    std::cout.precision(2);
    std::cout << std::fixed;
    std::cout << "Position: " << x
              << "\t\tFlux Value: " << y;
    
    std::cout.precision(5);
    std::cout << "\t\tFE Value: " << sum
              << "\t\tDifference: " << std::scientific << y - sum << std::endl;
}
double time = (double)(clock() - tStart)/CLOCKS_PER_SEC;
myfile << "Number of Particle: " + std::to_string(basis.N) + "\n";
myfile << "Time take: " + std::to_string(time) + "\n";
for (int m=0; m < basis.M; m++)
{
myfile << std::to_string(tally.coefficient_matrix[m]) + " * P_" + std::to_string(m) + "(x) +/- " + std::to_string(tally.unc_matrix[m]) + "With an R^2 value of " + std::to_string(tally.R_sqr_value[m]) + "\n";
}
myfile << "\n";
myfile.close();
}
