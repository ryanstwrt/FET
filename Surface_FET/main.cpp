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

int main (int argc, char** argv)
{
ofstream myfile;
ifstream input;
input.open("/opt/Shift_inputs/point10e3.rst");
myfile.open ("misc.txt", ios::in | ios::app);
clock_t tStart = clock ();

const std::size_t poly_order = 10;
const std::size_t poly_terms = poly_order + 1;
const std::size_t num_surfaces = 1;
double xs;
double wt;
double x1;
double y1;
double z1;
double counter;
double counter1;
double counter2;
char c;
std::string::size_type sz;
std::string str;

tally_info tally (poly_order);
legendre_info basis(poly_order, num_surfaces);
particle_info a;
basis.surface_index=0;
FET_solver solver;

while( !input.eof() )
{
	getline( input, str);
	int size_old = a.size;
	a.size = str.size();	
	if(a.size <= 3 && size_old >= 3)
	{
		a.b_alive = 0;
		solver.collision_eval2 (basis, a, poly_terms);
	}
	else if(a.size <= 12 && a.size > 3)
	{
		a.b_alive = 1;
		a.xs_tot = std::stod(str,&sz);
	}
	else if(a.size > 12)
	{
		a.b_alive = 1;
		std::istringstream iss(str);
		std::vector<std::string> temp((std::istream_iterator<std::string>(iss)), std::istream_iterator<std::string>());
		a.b_weight = std::stod(temp[0],&sz); 
		a.x = std::stod(temp[1],&sz); 
		a.y = std::stod(temp[2],&sz); 
		a.z = std::stod(temp[3],&sz);

		solver.collision_eval (basis, a, poly_terms);
	}

}
std::cout<< basis.n_counter[0] << std::endl;
solver.get_current(basis, tally, poly_terms);

for (int m=0; m < poly_terms; ++m)
{
   for(int n=0; n<poly_terms; ++n)
    {
      for(int i=0; i<poly_terms; ++i)
      {
	if(/*tally.flux_unc_matrix[m][n][i] >= tally.flux_unc_matrix[m][n][i]/100 || */tally.flux_R_matrix[m][n][i] >= 1)
	{
	std::cout<<"P("<<m<<")("<<n<<")("<<i<<") = " << tally.flux_matrix[m][n][i] << " +/- " << tally.flux_unc_matrix[m][n][i] << " w/ " << tally.flux_R_matrix[m][n][i] <<std::endl;
	if(tally.flux_R_matrix[m][n][i] >= 10)
	{
	counter++;
	}
	counter1++;
	}
//	if(tally.flux_R_matrix[m][n][i] > 1)
	counter2++;
       }
	std::cout<<std::endl;
    }
}
std::cout<<"Percentage of R^2 > 1: " <<counter/counter2<<std::endl;
std::cout<<"Percentage of R^2 > 10: " <<counter1/counter2<<std::endl;

double time = (double)(clock() - tStart)/CLOCKS_PER_SEC;
myfile << "Number of Particle: " + std::to_string(basis.n_counter[0]) + "\n";
myfile << "Time take: " + std::to_string(time) + "\n";
myfile << "Percentage of R^2 > 1: " + std::to_string(counter/counter2) +"\n";
myfile << "Percentage of R^2 > 10: " + std::to_string(counter1/counter2) +"\n";

for (int m=0; m < poly_terms; ++m)
{
   for(int n=0; n<poly_terms; ++n)
    {
      for(int i=0; i<poly_terms; ++i)
      {
	myfile << "+";
	myfile << tally.flux_matrix[m][n][i];
	myfile << "*LegendreP[";
	myfile << std::to_string(m);
	myfile << ",x]*LegendreP[";
	myfile << std::to_string(n);
	myfile << ",y]*LegendreP[";
	myfile << std::to_string(i);
	myfile << ",z]";
       }
    }
}

myfile << "\n";
myfile << "\n";
myfile.close();
input.close();

}

