#include<iostream>
#include<stdlib.h>

using namespace std;

int main ()
{

system("git add FET_legendre_solver.cc FET_post_process.cc FET.cc FET.hh main.cpp time.txt ../Report/Report.tex");
cout<<"Please write the message for update: "<<endl;
char message_temp [200];
cin.getline(message_temp,200);

string message =  "git commit -m \" ";
message += message_temp;
message += "\" ";
system(message.c_str());
system("git push origin master");

}
