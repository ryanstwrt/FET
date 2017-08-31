
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
#include <iterator>
#include <sstream>
using namespace  std;

#include"FET.hh"
#include "Distribution.hh"

int main (int argc, char** argv)
{
ofstream myfile;
ifstream input;
//input.open("/opt/Shift_inputs/shift_input_cell_small.rst");
myfile.open ("time.txt", ios::in | ios::app);
clock_t tStart = clock ();

const std::size_t poly_order = 5;
const std::size_t poly_terms = poly_order + 1;
const std::size_t N = 1e6;
const std::size_t num_surfaces = 1;
double xs;
double wt;
double x1;
double y1;
double z1;
char c;
std::string::size_type sz;
std::string str;
Distribution fluxshape;

tally_info tally (poly_order);
legendre_info basis(poly_order, num_surfaces);
particle_info a;
basis.surface_index=0;
FET_solver solver;
/*
while( !input.eof() )
{
getline( input, str);
int size = str.size();
//std::cout<<str<<std::endl;
if(size <= 3)
{
a.b_alive = 0;
}
else if(size <= 12 && size > 3)
{

a.xs_tot = std::stod(str,&sz);
}
else if(size > 12)
{
std::istringstream iss(str);
std::vector<std::string> temp((std::istream_iterator<std::string>(iss)), std::istream_iterator<std::string>());

	a.b_weight = std::stod(temp[0],&sz); //std::stod(str,&sz);
	a.x = std::stod(temp[1],&sz); //std::stod (str.substr(sz));
	a.y = std::stod(temp[2],&sz); //std::stod (str.substr(str.substr(sz)));
	a.z = std::stod(temp[3],&sz); //std::stod (str.substr(sz));
}

	solver.collision_eval (basis, a, poly_terms);
	solver.surface_eval (basis, a, poly_terms);

}
*/

for(int n = 0; n<N; n++)
{
	a.get_particle(basis, a);

	solver.collision_eval (basis, a, poly_terms);
	solver.collision_eval2 (basis, a, poly_terms);
}

solver.get_current(basis, tally, poly_terms, N);
/*
for(int s=0; s<basis.num_surfaces; s++)
{	
	std::cout<<"Surface Number "<<s<<std::endl;
	for(int m=0; m<poly_terms; m++)
	{
		std::cout<<tally.current_matrix_x[m]<<" * P_" <<m<<"(x) +/- "<<tally.current_unc_matrix[m]<<"    ";
		std::cout<<"With an R^2 value of "<<tally.current_R_sqr_value[m]<<std::endl;
	}

std::cout<<std::endl;
}*/

for(int s=0; s<basis.num_surfaces; s++)
{	
	std::cout<<"Surface Number "<<s<<std::endl;
	for(int m=0; m<poly_terms; m++)
	{
		std::cout<<tally.flux_matrix_x[m]<<" * P_" <<m<<"(x) +/- "<<tally.flux_unc_matrix[m]<<"    ";
		std::cout<<"With an R^2 value of "<<tally.flux_R_sqr_value[m]<<std::endl;
	}

std::cout<<std::endl;
}
/*
for(int cp = 0; cp <= 10; ++cp)
{
 double x = double(cp) / 10 * 2 -1;
 double y = double(cp) / 10 * 2 -1;
 double z = fluxshape(x) * fluxshape(y);
 double sum = 0;
 double sum1 = 0;
 double sum_tot;
 basis.Pn(poly_terms, x);

  for(int m = 0; m < poly_terms; m++)
  {
	sum += tally.current_matrix_x[m] * basis.P_n[m];
	sum1 += tally.current_matrix_y[m] * basis.P_n[m];
	sum_tot = sum*sum1;
  }

    std::cout.precision(2);
    std::cout << std::fixed;
    std::cout << "Position x: " << x
	      << "\t Position y: " << y
              << "\t Flux Value: " << z;
    
    std::cout.precision(3);
    std::cout << std::fixed;
    std::cout << "\t FE Value: " << sum_tot
              << "\t Difference: " << std::scientific << z - sum_tot <<std::endl;
}
*/

for(int cp = 0; cp <= 10; ++cp)
{
 double x = double(cp) / 10 * 2 -1;
 double y = double(cp) / 10 * 2 -1;
 double z = double(cp) / 10 * 2 -1;
 double t = fluxshape(x) * fluxshape(y) * fluxshape(z);
 double sum = 0;
 double sum1 = 0;
 double sum2 = 0;
 double sum_tot;
 basis.Pn(poly_terms, x);

  for(int m = 0; m < poly_terms; m++)
  {
	sum += tally.flux_matrix_x[m] * basis.P_n[m];
	sum1 += tally.flux_matrix_y[m] * basis.P_n[m];
	sum2 += tally.flux_matrix_z[m] * basis.P_n[m];
	sum_tot = sum*sum1*sum2;
  }

    std::cout.precision(2);
    std::cout << std::fixed;
    std::cout << "Position x: " << x
	      << "\t Position y: " << y
	      << "\t Position y: " << z
              << "\t Flux Value: " << t;
    
    std::cout.precision(3);
    std::cout << std::fixed;
    std::cout << "\t FE Value: " << sum_tot
              << "\t Difference: " << std::scientific << t - sum_tot <<std::endl;
}

double time = (double)(clock() - tStart)/CLOCKS_PER_SEC;
myfile << "Number of Particle: " + std::to_string(N) + "\n";
myfile << "Time take: " + std::to_string(time) + "\n";
for (int m=0; m < poly_terms; ++m)
{
myfile << std::to_string(tally.current_matrix_x[m]) + " * P_" + std::to_string(m) + "(x) +/- " + std::to_string(tally.current_unc_matrix[m]) + " With an R^2 value of " + std::to_string(tally.current_R_sqr_value[m]) + "\n";
}
myfile << "\n";
myfile.close();
}

