#include <iostream>
#include "fitsimage.h"


int main(int nargs, char**args)
{
  FitsImage iml(args[1]);
  FitsImage imr(args[2]);

  int nx = iml.Nx(); 
  int ny = iml.Ny(); 
  Image imtot(iml.Nx()+imr.Nx(),ny);


  for (int j=0; j< ny; ++j)
    for (int i=0;i<nx;++i)
      {
	imtot(i,j) = iml(i,j);
      }

  for (int j=0; j< ny; ++j)
    for (int i=0;i<nx;++i)
      {
	imtot(i+nx,j) = imr(i,j);
      }

  FitsImage out(args[3],iml,imtot);
}
