#ifndef PIXELBLOCK__H
#define PIXELBLOCK__H

#include <string>
#include <string.h>
#include "array2d.h"

typedef double PixelType; 


/* This class implements algebra on PixelBlocks which could be inside
   Array2d. They are here because they require both a Copy contructor
   and an operator =, which are hard to put in Array2D, or would be
   inefficient because we could then not make use of the fact that we
   deal with a "simple type", as we do here. There is perhaps a way of
   specializing a template for simple types, but I don't know it. If
   we need an other class with the same algebra capabilities, we can
   derive the Array2D template for simple types (Array2DSimple) , add
   into it the algebra functionnalities, and derive here from
   Array2DSimple. */

/*! a class that holds pixels of type PixelType (typically float),
  within arbitrary bounds (at variance with Image). Read and Write to
  Fits capabilities, plus algebra. */
class PixelBlock : public Array2D<PixelType>
{
  public :

    PixelBlock() : Array2D<PixelType>() {};

    PixelBlock(const IntFrame &I) : Array2D<PixelType>(I) {};

  PixelBlock(const int Imin, const int Jmin, const int Imax, const int Jmax):
  Array2D<PixelType>(Imin, Jmin, Imax, Jmax) {};

    /*! read pixels. The read stamp is defined via the contained IntFrame,
      with offsets provided through XRef and YRef. Assumes the PixelBlock
      has been Allocate'd.*/
    bool ReadFromFits(const std::string &FileName, 
		      const int XRef, const int YRef);


    bool WriteFits(const std::string &FileName) const;


   PixelBlock operator - (const PixelBlock &Right) const;

   void operator -= (const PixelBlock &Right);

   void operator -= (const double &Val)
   {
     const double* pend = this->Data() + Ntot();
     for (double* p = this->Data(); p < pend; ++p) (*p) -= Val;
   }     


   PixelBlock( const PixelBlock &Right)
   {
     *this = Right;
   }

   void operator = (const PixelBlock &Right)
   {
     Allocate((IntFrame &)(Right));
     memcpy(this->Data(), Right.Data(), Ntot()*sizeof(PixelType));
   }

   //! position of the pixel w.r.t beginning
  int PixelIndex(const int I, const int J) const
  {
#ifdef SIMPHOTFIT_CHECK_BOUNDS
    if (I<xmin || I >= xmax || J<ymin || J >= ymax)
      {
	std::cout << "index out of range in ModelPix::PixelIndex  (i,j,limits) : "
		  << I << ',' << J << ' ' << IntFrame(*this) << std::endl;
        abort();
      }
#endif
    return ((I-xmin)+(J-ymin)*nx);
  }

  void Moments(double &Xmean, double &Ymean, 
	       double &VarXX , double &VarYY, double &VarXY) const;


};


double ScalProd(const PixelBlock &, const PixelBlock &);


void GaussianFill(PixelBlock &K, const double SigX, const double SigY, const double dx=0, const double dy=0);



#endif /* PIXELBLOCK__H */
