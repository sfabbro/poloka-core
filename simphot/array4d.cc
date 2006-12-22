#include "array4d.h"


#include <algorithm>

#define BIGINT 1000000

size_t Array4D::MemSize() const
{
  size_t tot = 0;
  for (int y = ymin; y <ymax; ++y)
    for (int x = xmin; x < xmax; ++x)  
      {
	tot += (*this)(x,y).Ntot();
      }
  tot *= sizeof(CoeffType); // words to bytes
  tot += sizeof(CoeffBlock)*Ntot(); // add administrative data
  return tot;
}
  

using namespace std; // for min and max

//! Compute the boundary that contains all subblocks
IntFrame Array4D::InternalFrame() const
{
  int imin = BIGINT;
  int imax = -imin;
  int jmin = imin;
  int jmax= -jmin;
  for (int y = ymin; y <ymax; ++y)
    for (int x = xmin; x < xmax; ++x)
      {
	const CoeffBlock& block = (*this)(x,y);
	if (block.Ntot() == 0) continue;
	imin = min(block.xmin, imin);
	jmin = min(block.ymin, jmin);
	imax = max(block.xmax, imax);
	jmax = max(block.ymax, jmax);
      }
  return IntFrame(imin,jmin,imax,jmax);
}


void Array4D::Transpose(Array4D &Result) const
{
/* Checked that Result.Transposed() == *(this) on non-trivial contents */
  int imin = BIGINT;
  int imax = -imin;
  int jmin = imin;
  int jmax= -jmin;
  for (int y = ymin; y <ymax; ++y)
    for (int x = xmin; x < xmax; ++x)
      {
	const CoeffBlock& block = (*this)(x,y);
	if (block.Ntot() == 0) continue;
	imin = min(block.xmin, imin);
	jmin = min(block.ymin, jmin);
	imax = max(block.xmax, imax);
	jmax = max(block.ymax, jmax);
      }
  Result.Allocate(imin,jmin,imax,jmax);

  // find extrema for each i,j

  for (int j = jmin; j <jmax; ++j)
  for (int i = imin; i <imax; ++i)
    {
      int amin = BIGINT;
      int amax = -BIGINT;
      int bmin = BIGINT;
      int bmax = -bmin;
      int nhits = 0;
      for (int b = ymin; b<ymax; ++b)
	for (int a = xmin; a<xmax; ++a)
	  {
	    const CoeffBlock &block = (*this)(a,b);
	    if (block.Ntot() == 0) continue;
	    if (block.IsInside(i,j))
	      {
		amin = min(a,amin);
		amax = max(a,amax);
		bmin = min(b,bmin);
		bmax = max(b,bmax);
		nhits++;
	      }
	  }
      if (nhits) 
	{
	  Result(i,j).Allocate(amin, bmin, amax + 1, bmax + 1);
	  Result(i,j).SetVal(0.);
	}
      // DEBUG
      // I have checked here that when not oversmapling, the blocks have the right size.
      // if (Result(i,j).Ntot() > 9) abort();
    }
  for (int b = ymin; b<ymax; ++b)
    for (int a = xmin; a<xmax; ++a)
	  {
  	    const CoeffBlock &block = (*this)(a,b);
	    for (int j = block.ymin; j < block.ymax; ++j)
	    for (int i = block.xmin; i < block.xmax; ++i)
	      {
		Result(i,j)(a,b) = block.at(i,j);
	      }
	  }	
}


void Array4D::Convolve(const Array2D<double> &Kernel, 
		       Array4D &Result) const
{
  Result.Allocate((const IntFrame &)(*this));
  int kxmin = Kernel.xmin;
  int kymin = Kernel.ymin;
  int kxmax = Kernel.xmax;
  int kymax = Kernel.ymax;

  for (int j = ymin; j < ymax; ++j)
    for (int i = xmin; i < xmax; ++i)
      {
	const CoeffBlock &bin = (*this)(i,j);
	CoeffBlock &bout = Result(i,j);
	// TODO incorporate integration limits that should be provided
	// as arguments (can be deduced from the "InternalFrame", but it 
	// takes time) . outFrame is too large!
	IntFrame outFrame(bin.xmin+Kernel.xmin, bin.ymin+Kernel.ymin,
			  bin.xmax+Kernel.xmax-1, bin.ymax+Kernel.ymax-1);
	bout.Allocate(outFrame);
	for (int d=bout.ymin; d < bout.ymax; ++d)
	  {
	    int bmin = max(d-kymax+1, bin.ymin);
	    int bmax = min(d-kymin+1, bin.ymax);
	    for (int c=bout.xmin; c < bout. xmax; ++c)
	      {
		double sum = 0;
		int amin = max(c-kxmax+1, bin.xmin);
		int amax = min(c-kxmin+1, bin.xmax);
#define FASTWAY
#ifdef FASTWAY
		for (int b = bmin; b<bmax;++b)
		  {
		    const CoeffType *pbin = &bin(amin,b);
		    const double *pk = &Kernel(c-amin, d-b);
		    for (int a = amin; a<amax;++a)
		      {
			sum += (*pbin)*(*pk);
			++pbin; --pk;
		      }
		  }
#else
		for (int b = bmin; b<bmax;++b)
		for (int a = amin; a<amax;++a)
		  sum += bin(a,b)*Kernel(c-a,d-b);
#endif
		bout(c,d) = sum;
	      }
	  }
      }
}

const CoeffBlock& Array4D::GetElement(const int i, const int j) const 
{ return (*this)(i,j);}
