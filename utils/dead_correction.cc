#include <iostream>
#include <string>

#include "fitsimage.h"

int main(int argc,char **args)
{
if (argc <=1)
  {
    cout << " dead_correction <fitsfile(s)> " << endl;
    exit(-1);
  }
  for (int i=1; i< argc; ++i)
    {
      string inname = args[i];
      FitsImage in(args[i]);
      if (!in.IsValid())
	{
	  cerr << " dead_correction : invalid file : "  << args[i] << endl;
	  continue;
	}
      string outname = "dc_" + inname;
      FitsImage out(outname,in,in);
      int nx = out.Nx();
      int ny = out.Ny();

      int j = 0;
      while (j<ny) 
	{
	  int i = 0;
	  while (i<nx)
	  {     
	    double val = ceil(0.25*(out(i,j) + out(i+1,j) + out(i,j+1) + out(i+1,j+1)));
	    out(i,j) = val;
	    out(i+1,j) = val;
	    out(i,j+1) = val;
	    out(i+1,j+1) = val;
	    i = i + 2;
	  }
	  j = j + 2;
	}
    }
  return 1;

}
