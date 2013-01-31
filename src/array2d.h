#ifndef ARRAY2D__H
#define ARRAY2D__H


#include <iostream>

#include "intframe.h"
#include <stdlib.h>

#define PIXEL_LOOP(block, I, J) \
  for (int J = block.ymin; J < block.ymax; ++J) \
  for (int I = block.xmin; I < block.xmax; ++I)




#define ARRAY2D_CHECK_BOUNDS

//! A container with contents indexed by 2 integers: a "generalized" Image
/*! it is generalized because contents may be anything and ranges do not
necessarily start at 0,0. The limits are handled via IntFrame. The convention
is that pixel indices start at xmin,ymin and end at xmax-1, ymax-1. */ 
template <class T> class Array2D : public IntFrame
{

 public: // the friend declarations are just useless.
  // we refer to coordinates using x and y because in the implementation
  // of Array4D, we use 2 sets of variables: a,b, and i,j
  // I hope it will eventually make the code almost understandable

 private:
  T* data; 
  T* data00;

  /* if we switch to a more elaborated allocator, we have to modify where
     ALLOCATOR is written. In the present state, copies are forbidden */
 public:

  Array2D() : 
    IntFrame()
  {
    data = NULL; // ALLOCATOR
    Allocate();
  }



  Array2D(int Xmin, int Ymin, int Xmax, int Ymax) : 
    IntFrame(Xmin,Ymin, Xmax, Ymax)
  {
    data = NULL; // ALLOCATOR
    Allocate();
  }

  Array2D(const IntFrame &Frame) : IntFrame(Frame)
  {
    data = NULL; // ALLOCATOR
    Allocate();
  }    

  void Allocate(const IntFrame &IF)
  {
    (IntFrame &)(*this)= IF;
    Allocate();
  }

  void Allocate(int Xmin,  int Ymin, int Xmax, int Ymax)
  {
    xmin = Xmin; xmax = Xmax; ymin =  Ymin;  ymax = Ymax;
    nx = xmax -xmin;
    Allocate();
  }

  void Allocate()
  {
    if (data) delete[] data; // ALLOCATOR
    data = data00 = NULL;
    int ntot = Ntot();
    if (ntot > 0)
      {
	data = new T[ntot]; // ALLOCATOR
	data00 = &data[(0-xmin)+(0-ymin)*nx];
      }
  }



  ~Array2D() // ALLOCATOR. may become useless if handled with a container
  {
    if (data) 
      delete [] data; 
  }

  T* Data() {return data;}
  const T* Data() const {return data;}
  T* begin() {return data;}
  const T* begin() const {return data;}
  T* end() { return data+Ntot();}
  const T* end() const { return data+Ntot();}


  inline const T& operator() (const int x, const int y) const
  {
#ifdef  ARRAY2D_CHECK_BOUNDS
    if (IsOutside(x,y))
      {
	std::cout << "index out of range in Array2D (i,j, interval)" 
		  << x << ',' << y << ',' << IntFrame(*this) << std::endl;
	abort();
      }
#endif
    return  data00[x+y*nx];
  }


  inline T& operator() (const int x, const int y)
  {
#ifdef  ARRAY2D_CHECK_BOUNDS
    if (x<xmin || x>= xmax || y < ymin || y>=ymax)
      {
	std::cout << "index out of range in Array2D (i,j, interval)" 
		  << x << ',' << y << ',' << IntFrame(*this) << std::endl;
	abort();
      }
#endif
    return  data00[x+y*nx];
  }

  //access without bound checking
  inline const T& at(const int x, const int y) const
  {
    return data00[x+y*nx];
  }

  inline T& at(const int x, const int y)
  {
    return data00[x+y*nx];
  }


  bool operator == (const Array2D<T> & Right) const
  {
    if (xmin != Right.xmin || ymin != Right.ymin 
	|| xmax != Right.xmax || ymax != Right.ymax) return false;
    for (int y =ymin; y < ymax; ++y)
      for (int x= xmin; x < xmax; ++x)
	if (!((*this)(x,y) == Right(x,y))) return false;
    return true;
  }

   
  void operator *= (const double &Fact)
  {
    const T* pend = data + Ntot();
    for (T* p = data; p < pend; ++p) (*p) *= Fact;
  }


  void operator -= (const double &Val)
  {
    const T* pend = data + Ntot();
    for (T* p = data; p < pend; ++p) (*p) -= Val;
  }

  void operator += (const double &Val)
  {
    const T* pend = data + Ntot();
    for (T* p = data; p < pend; ++p) (*p) += Val;
  }


  void SetVal(const T& Val)
  {
    T* pend = data + Ntot();
    for (T* p = data; p < pend; ++p) (*p) = Val;
  }


  T Sum() const
  {
    T sum = 0;
    T* pend = data + Ntot();
    for (T* p = data; p < pend; ++p) sum += *p;
    return sum;
  }  

 private:
  // no copy handling yet
  Array2D<T>(const Array2D<T> & A); // would have to use copy contructor of elements
  void operator = (const Array2D<T> & A); // would have to use operator = of elements

};


#endif /* ARRAY2D__H */
