// -*- C++ -*-
#include <iostream>
#include <iomanip>
#include <string>

int main(int argc, char** argv)
{
  string f1, f2, f3;
  f1="pomme"; f2="poire"; f3="abricot";
  
  //  std::cout.flags(ios::left);
  //  std::cout.width(10);
  
  std::cout << "|" << setw(10) << f1
	    << "|" << setw(10) << f2
	    << "|" << setw(10) << f3
	    << "|"
	    << endl;
  
}

