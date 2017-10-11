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
input.open("/opt/Shift_inputs/fuel_rod1.5.rst");
myfile.open ("fuel_rod.txt", ios::in | ios::app);
clock_t tStart = clock ();

FET_solver solver;
initial_info info;
solver.initializer(info);
tally_info tally (info.poly_order, info.num_tallies);
legendre_info basis (info.poly_order, info.num_tallies);
particle_info a;

double xs;
double wt;
double x1;
double y1;
double z1;
char c;
int counter = 0;
std::string::size_type sz;
std::string str;

while( !input.eof() )
{
	counter++;
	getline( input, str);
	int size_old = a.size;
	a.size = str.size();
	
	if(a.size <= 3 && size_old >= 3)
	{
		a.alive = 0;
		solver.collision_eval2 (basis, a, info.poly_terms);
	}
	else if(a.size <= 12 && a.size > 3)
	{
		a.alive = 1;
		a.xs_tot = std::stod(str,&sz);
	}
	else if(a.size > 12)
	{
		a.alive = 1;
		std::istringstream iss(str);
		std::vector<std::string> temp((std::istream_iterator<std::string>(iss)), std::istream_iterator<std::string>());
		a.wt = std::stod(temp[0],&sz); 
		a.x = std::stod(temp[1],&sz); 
		a.y = std::stod(temp[2],&sz); 
		a.z = std::stod(temp[3],&sz);

		solver.collision_eval (basis, a, info.poly_terms);
	}
}

std::cout<< basis.n_counter[0] << std::endl;
solver.get_current(basis, tally, info.num_tallies, info.poly_terms);
solver.cleanup(tally, info.poly_terms, info.num_tallies);

//To Do: Remove coefficients with greater than 10, and warn the user for coefficients greater than 1


for(int k=0; k<info.num_tallies; ++k)
{
//std::cout<<tally.R_greater[k]<<"  "<<tally.R_greater[k+1]<<"  "<<std::endl;
std::cout<<"Percentage of R^2 > 1: " <<tally.R_greater[2*k]/tally.total_coeff[k]*100<<std::endl;
std::cout<<"Percentage of R^2 > 10: " <<tally.R_greater[2*k+1]/tally.total_coeff[k]*100<<std::endl;
std::cout<<basis.n_counter[k]<<std::endl;
std::cout<<std::endl;
}

double time = (double)(clock() - tStart)/CLOCKS_PER_SEC;
myfile << "Time take: " + std::to_string(time) + "\n";

for(int k=0; k<info.num_tallies; ++k)
{
if(k==0)
myfile<<"Fuel Rod Only \n";

if(k==1)
myfile<<"Fuel + Water";

  myfile << "Tally " << k << " \n";
  myfile << "Number of Particle: " + std::to_string(basis.n_counter[k]) + "\n";
  myfile << "Percentage of R^2 > 1: " + std::to_string((tally.R_greater[2*k])/tally.total_coeff[k]*100) +"\n";
  myfile << "Percentage of R^2 > 10: " + std::to_string(tally.R_greater[2*k+1]/tally.total_coeff[k]*100) +"\n";
  for (int m=0; m < info.poly_terms; ++m)
  {
    for(int n=0; n<info.poly_terms; ++n)
    {
      for(int i=0; i<info.poly_terms; ++i)
      {

	myfile << "+";
	myfile << tally.flux_matrix[m][n][i][k];
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
}

myfile << "\n";
myfile << "\n";
myfile.close();
input.close();

}

