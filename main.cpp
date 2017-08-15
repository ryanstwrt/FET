
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

const std::size_t poly_order = 5;
const std::size_t poly_terms = poly_order + 1;
const std::size_t N = 1e6;
const std::size_t num_surfaces = 1;

Distribution fluxshape;

tally_info tally;
legendre_info basis(poly_order, num_surfaces);
particle_info a;
tally.surface_index=0;
a.a_n.resize(poly_terms, 0.0);

initialize_tally_info (tally, poly_terms);

for(int n = 0; n<N; n++)
{
	a.get_particle(basis, a);

	surface_eval (basis, a, poly_terms);
}

get_current(basis, tally, poly_terms, N);

for(int s=0; s<tally.num_surfaces; s++)
{	
	std::cout<<"Surface Number "<<s<<std::endl;
	for(int m=0; m<poly_terms; m++)
	{
		std::cout<<tally.current_matrix[m]<<" * P_" <<m<<"(x) +/- "<<tally.unc_matrix[m]<<"    ";
		std::cout<<"With an R^2 value of "<<tally.R_sqr_value[m]<<std::endl;
	}

std::cout<<std::endl;
}

for(int cp = 0; cp <= 10; ++cp)
{
 double x = double(cp) / 10 * 2 -1;
 double y = fluxshape(x);
 double sum = 0;
 basis.Pn(poly_terms, x);

  for(int m = 0; m < poly_terms; m++)
	sum += tally.current_matrix[m] * basis.P_n[m];
    std::cout.precision(2);
    std::cout << std::fixed;
    std::cout << "Position: " << x
              << "\t\tFlux Value: " << y;
    
    std::cout.precision(5);
    std::cout << "\t\tFE Value: " << sum
              << "\t\tDifference: " << std::scientific << y - sum << std::endl;
}
double time = (double)(clock() - tStart)/CLOCKS_PER_SEC;
myfile << "Number of Particle: " + std::to_string(N) + "\n";
myfile << "Time take: " + std::to_string(time) + "\n";
for (int m=0; m < poly_terms; ++m)
{
myfile << std::to_string(tally.current_matrix[m]) + " * P_" + std::to_string(m) + "(x) +/- " + std::to_string(tally.unc_matrix[m]) + " With an R^2 value of " + std::to_string(tally.R_sqr_value[m]) + "\n";
}
myfile << "\n";
myfile.close();
}
