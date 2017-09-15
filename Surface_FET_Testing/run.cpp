#include<iostream>
#include<stdlib.h>

using namespace std;

int main ()
{

system("g++ FET_legendre_solver.cc FET_post_process.cc FET.cc FET.hh main.cpp Distribution.cpp Distribution.hh -o test -std=c++11");
system("./test");

}
