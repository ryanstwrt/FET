

#ifndef FET_H
#define FET_H
#include"FET.hh"
#endif

int main ()
{

legendre_info basis1;
legendre_info basis2;
tally_info first;
particle_info a;
particle_info b;

std::cout<<"test"<<std::endl;

for(int surf=0; surf<2;surf++)
{
if (surf == 0)
{
std::cout<<"test"<<std::endl;
surface_eval(basis1,a);
std::cout<<"test"<<std::endl;
}
else
{
surface_eval(basis2,b);
}
std::cout<<"The scaled current for surface "<<surf<<" is described as: "<<std::endl;
for(int n=0; n<basis1.M; n++)
{	
if (surf == 0)
{
	std::cout<<basis1.current[n]<<"*P_"<<n<<"(x) +/- "<<basis1.var_a_n[n]<<std::endl;
}
else
{
	std::cout<<basis2.current[n]<<"*P_"<<n<<"(x) +/- "<<basis2.var_a_n[n]<<std::endl;

}
}

std::cout<<std::endl;
}
}
