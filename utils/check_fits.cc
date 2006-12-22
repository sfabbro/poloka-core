#include <iostream>
#include "fitsimage.h"
#include "fitsexception.h"



int main(int nargs, char **args)
{
  int status;
  for (int i=1; i< nargs; ++i)
    {
      try 
	{
	  FitsImage toto(args[i]);
	}
      catch (FitsException toto)
	{
	  std::cout << " BAD_FILE : " << args[i] << std::endl;
	  status = 1;
	}
    }
  return status;
}
