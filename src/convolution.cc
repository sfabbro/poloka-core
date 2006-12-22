#include "convolution.h"
#include "exceptions.h"
#include <cmath> // for exp()

#ifndef M_PI
#define     M_PI            3.14159265358979323846  /* pi */
#endif


Image Convole_MakeNoyau(double r, int l)
{ 
  int d = 2*l+1;
  Image noyau(d,d);
  TRY{
  {for (int i=0; i<d; i++)
    {for (int j=0; j<d; j++)
      noyau(i,j) = (float) (exp( -( (i - l)*(i -l) + (j- l)*(j -l) )/ (2 * r * r ) ) );
    }}
   }
 CATCHALL{THROW_SAME;}
 ENDTRY
  return noyau;
}

double * Convole_MakeNoyau_1D(double r, int l)
{ 
  double * p; 
  int d = 2*l+1;
  TRY{
  p = new double[d] ;
  for (int i=0; i<d; i++)
    {
      double x = ((double) (i-l)) ;
      p[i] = exp( -0.5 * (x * x)/(r*r));
      /*double u1 = (x-0.5)/(sqrt(2.)*r) ;
      double u2 = (x+0.5)/(sqrt(2.)*r) ;
      p[i] = erf(u2)-erf(u1) ;*/

    }
  }
  CATCHALL{THROW_SAME;}
  ENDTRY
    /*for (int i=0;i <d;++i)
      { 
      double x = ((double) (i-l)) ;
      cout <<  exp( -0.5 * (x * x)/(r*r))  << " " << p[i]/p[l] << endl  ;
      }*/
  return p;
}



Image Convole_MakeNoyau(double rx, double ry , int l)
{ 
  int d = 2*l+1;
  Image noyau(d,d);
  TRY{
  {for (int y=0; y<d; y++)
       {for (int x=0; x<d; x++)
            noyau(x,y) = (float)   (exp( -(  (x - l)*(x -l) / (2 * rx * rx ) + 
					     (y- l)*(y -l) / (2 * ry * ry )  ) ) );
       }
  }
   }
  CATCHALL{THROW_SAME;}
  ENDTRY
    return noyau;
  } 


#define Convole_LDEMICOTE Convole_2 
#define Convole1DX_LDEMICOTE Convole1DX_2 
#define Convole1DY_LDEMICOTE Convole1DY_2
#define LDEMICOTE 2
#include "convtemplate.h"
#undef  Convole_LDEMICOTE
#undef Convole1DX_LDEMICOTE 
#undef Convole1DY_LDEMICOTE 
#undef LDEMICOTE

#define Convole_LDEMICOTE Convole_3 
#define Convole1DX_LDEMICOTE Convole1DX_3 
#define Convole1DY_LDEMICOTE Convole1DY_3
#define LDEMICOTE 3
#include "convtemplate.h"
#undef  Convole_LDEMICOTE
#undef Convole1DX_LDEMICOTE 
#undef Convole1DY_LDEMICOTE 
#undef LDEMICOTE

#define Convole_LDEMICOTE Convole_4 
#define Convole1DX_LDEMICOTE Convole1DX_4
#define Convole1DY_LDEMICOTE Convole1DY_4
#define LDEMICOTE 4
#include "convtemplate.h"
#undef  Convole_LDEMICOTE
#undef Convole1DX_LDEMICOTE 
#undef Convole1DY_LDEMICOTE 
#undef LDEMICOTE


#define Convole_LDEMICOTE Convole_5 
#define Convole1DX_LDEMICOTE Convole1DX_5
#define Convole1DY_LDEMICOTE Convole1DY_5
#define LDEMICOTE 5
#include "convtemplate.h"
#undef  Convole_LDEMICOTE
#undef Convole1DX_LDEMICOTE 
#undef Convole1DY_LDEMICOTE 
#undef LDEMICOTE


#define Convole_LDEMICOTE Convole_6 
#define Convole1DX_LDEMICOTE Convole1DX_6
#define Convole1DY_LDEMICOTE Convole1DY_6
#define LDEMICOTE 6
#include "convtemplate.h"
#undef  Convole_LDEMICOTE
#undef Convole1DX_LDEMICOTE 
#undef Convole1DY_LDEMICOTE 
#undef LDEMICOTE

#define Convole_LDEMICOTE Convole_7 
#define Convole1DX_LDEMICOTE Convole1DX_7
#define Convole1DY_LDEMICOTE Convole1DY_7
#define LDEMICOTE 7
#include "convtemplate.h"
#undef  Convole_LDEMICOTE
#undef Convole1DX_LDEMICOTE 
#undef Convole1DY_LDEMICOTE 
#undef LDEMICOTE

#define Convole_LDEMICOTE Convole_8 
#define Convole1DX_LDEMICOTE Convole1DX_8
#define Convole1DY_LDEMICOTE Convole1DY_8
#define LDEMICOTE 8
#include "convtemplate.h" 
#undef  Convole_LDEMICOTE
#undef Convole1DX_LDEMICOTE 
#undef Convole1DY_LDEMICOTE 
#undef LDEMICOTE

#define Convole_LDEMICOTE Convole_10
#define Convole1DX_LDEMICOTE Convole1DX_10
#define Convole1DY_LDEMICOTE Convole1DY_10
#define LDEMICOTE 10
#include "convtemplate.h" 
#undef  Convole_LDEMICOTE
#undef Convole1DX_LDEMICOTE 
#undef Convole1DY_LDEMICOTE 
#undef LDEMICOTE

#define Convole_LDEMICOTE Convole_11
#define Convole1DX_LDEMICOTE Convole1DX_11
#define Convole1DY_LDEMICOTE Convole1DY_11
#define LDEMICOTE 11
#include "convtemplate.h" 
#undef  Convole_LDEMICOTE
#undef Convole1DX_LDEMICOTE 
#undef Convole1DY_LDEMICOTE 
#undef LDEMICOTE

#define Convole_LDEMICOTE Convole_12
#define Convole1DX_LDEMICOTE Convole1DX_12
#define Convole1DY_LDEMICOTE Convole1DY_12
#define LDEMICOTE 12
#include "convtemplate.h" 
#undef  Convole_LDEMICOTE
#undef Convole1DX_LDEMICOTE 
#undef Convole1DY_LDEMICOTE 
#undef LDEMICOTE

#define Convole_LDEMICOTE Convole_13
#define Convole1DX_LDEMICOTE Convole1DX_13
#define Convole1DY_LDEMICOTE Convole1DY_13
#define LDEMICOTE 13
#include "convtemplate.h" 
#undef  Convole_LDEMICOTE
#undef Convole1DX_LDEMICOTE 
#undef Convole1DY_LDEMICOTE 
#undef LDEMICOTE

#define Convole_LDEMICOTE Convole_14
#define Convole1DX_LDEMICOTE Convole1DX_14
#define Convole1DY_LDEMICOTE Convole1DY_14
#define LDEMICOTE 14
#include "convtemplate.h" 
#undef  Convole_LDEMICOTE
#undef Convole1DX_LDEMICOTE 
#undef Convole1DY_LDEMICOTE 
#undef LDEMICOTE


#define Convole_LDEMICOTE Convole_15
#define Convole1DX_LDEMICOTE Convole1DX_15
#define Convole1DY_LDEMICOTE Convole1DY_15
#define LDEMICOTE 15
#include "convtemplate.h" 
#undef  Convole_LDEMICOTE
#undef Convole1DX_LDEMICOTE 
#undef Convole1DY_LDEMICOTE 
#undef LDEMICOTE

#define Convole_LDEMICOTE Convole_20
#define Convole1DX_LDEMICOTE Convole1DX_20
#define Convole1DY_LDEMICOTE Convole1DY_20
#define LDEMICOTE 20
#include "convtemplate.h" 
#undef  Convole_LDEMICOTE
#undef Convole1DX_LDEMICOTE 
#undef Convole1DY_LDEMICOTE 
#undef LDEMICOTE

int Demi_Fenetre(int l, double r, double precision)
{
 int n = l ;

 if ( l < 0)
   {
     n = Demi_Fenetre_Optimale(r,precision);
   }


 if (n < 2 )
   n = 2 ;
 if ( n > 15 && n < 20 )
   n = 20 ;  
 if ( n > 20 )
   {
     cerr << "fenetre trop grande " << n << ", taking 20 " << endl ;
     n = 20 ;  
   }
  return(n);
}
      

void Convole(Image const& src, double r, int l , 
 Image const& mask, Image & a )
{
  int n = Demi_Fenetre(l,r);
TRY{
 switch (n)
     {
      case 2 :
      Convole_2( src,  r, mask,  a );
      break;
      case 3 :
      Convole_3(src,  r, mask,  a );
      break;
      case 4 :
      Convole_4(src,  r, mask,  a );
      break;
      case 5 :
      Convole_5(src,  r, mask,  a );
      break;
      case 6 :
      Convole_6(src,  r, mask,  a );
      break;
      case 7 :
      Convole_7(src,  r, mask,  a );
      break;
      case 8 :
      Convole_8(src,  r, mask,  a );
      break;
      case 10 :
      Convole_10(src,  r, mask,  a );
      break;
      case 11 :
      Convole_11(src,  r, mask,  a );
      break;
      case 12 :
      Convole_12(src,  r, mask,  a );
      break;
      case 13 :
      Convole_13(src,  r, mask,  a );
      break;
      case 14 :
      Convole_14(src,  r, mask,  a );
      break;
      case 15 :
      Convole_15(src,  r, mask,  a );
      break;
      case 20 :
      Convole_20(src,  r, mask,  a );
      break;
      default :
      Convole_20(src,  r,   mask,  a );

     }
 }
  CATCHALL{THROW_SAME;}
 ENDTRY
}

void Convole(Image const& src, double r, int l , Image& a )
{
  int n = Demi_Fenetre(l,r);

  TRY{
  switch (n)
     {
      case 2 :
      Convole_2( src,  r,   a );
      break;
      case 3 :
      Convole_3( src,  r,   a );
      break;
      case 4 :
      Convole_4( src,  r,  a );
      break;
      case 5 :
      Convole_5(src,  r,  a );
      break;
      case 6 :
      Convole_6(src,  r,  a );
      break;
      case 7 :
      Convole_7(src,  r,  a );
      break;
      case 8 :
      Convole_8(src,  r,  a );
      break;
      case 10 :
      Convole_10(src,  r,  a );
      break;
      case 11 :
      Convole_11(src,  r,  a );
      break;
      case 12 :
      Convole_12(src,  r,  a );
      break;
      case 13 :
      Convole_13(src,  r,  a );
      break;
      case 14 :
      Convole_14(src,  r,  a );
      break;
      case 15 :
      Convole_15(src,  r,  a );
      break;
      case 20 :
      Convole_20(src,  r,  a );
      break;
      default :
      Convole_20(src,  r,  a );

     }
  }
    CATCHALL{THROW_SAME;}
  ENDTRY

}

void Convole(Image const& src, double rx,double ry, int l , 
	     Image const& mask, Image& a )
{
  double r = max(rx,ry) ;
  int n = Demi_Fenetre(l,r);


 TRY{
 switch (n)
     {
      case 2 :
      Convole_2( src,  rx , ry ,  mask,  a );
      break;
      case 3 :
      Convole_3(src,  rx , ry ,  mask,  a );
      break;
      case 4 :
      Convole_4(src,  rx , ry ,  mask,  a );
      break;
      case 5 :
      Convole_5(src,  rx , ry ,  mask,  a );
      break;
      case 6 :
      Convole_6(src,  rx , ry ,  mask,  a );
      break;
      case 7 :
      Convole_7(src,  rx , ry ,  mask,  a );
      break;
      case 8 :
      Convole_8(src,  rx , ry ,  mask,  a );
      break;
      case 10 :
      Convole_10(src,  rx , ry ,  mask,  a );
      break;
      case 11 :
      Convole_11(src,  rx , ry ,  mask,  a );
      break;
      case 12 :
      Convole_12(src,  rx , ry ,  mask,  a );
      break;
      case 13 :
      Convole_13(src,  rx , ry ,  mask,  a );
      break;
      case 14 :
      Convole_14(src,  rx , ry ,  mask,  a );
      break;
      case 15 :
      Convole_15(src,  rx , ry ,  mask,  a );
      break;
      case 20 :
      Convole_20(src,  rx , ry ,  mask,   a );
      break;
      default :
      Convole_20(src,  rx , ry ,    mask,  a );

     }
 }
 CATCHALL{THROW_SAME;}
 ENDTRY

}

void Convole(Image const& src, double rx , double ry ,  int l , Image& a )
{
  double r = max(rx,ry) ;
  int n = Demi_Fenetre(l,r);

TRY{
  switch (n)
     {
      case 2 :
      Convole_2( src,  rx , ry ,    a );
      break;
      case 3 :
      Convole_3( src,  rx , ry ,    a );
      break;
      case 4 :
      Convole_4( src,  rx , ry ,   a );
      break;
      case 5 :
      Convole_5(src,  rx , ry ,   a );
      break;
      case 6 :
      Convole_6(src,  rx , ry ,   a );
      break;
      case 7 :
      Convole_7(src,  rx , ry ,   a );
      break;
      case 8 :
      Convole_8(src,  rx , ry ,   a );
      break;
      case 10 :
      Convole_10(src,  rx , ry ,   a );
      break;
      case 11 :
      Convole_11(src,  rx , ry ,   a );
      break;
      case 12 :
      Convole_12(src,  rx , ry ,   a );
      break;
      case 13 :
      Convole_13(src,  rx , ry ,   a );
      break;
      case 14 :
      Convole_14(src,  rx , ry ,   a );
      break;
      case 15 :
      Convole_15(src,  rx , ry ,   a );
      break;
      case 20 :
      Convole_20(src,  rx , ry ,   a );
      break;
      default :
      Convole_20(src,  rx , ry ,   a );

     }
  }
  CATCHALL{THROW_SAME;}
  ENDTRY

}

void Convole1DX (Image const& src, double r,  int l , Image& a )
{
  
  int n = Demi_Fenetre(l,r);
 TRY{
  switch (n)
     {
      case 2 :
      Convole1DX_2( src,  r ,    a );
      break;
      case 3 :
      Convole1DX_3( src,  r,    a );
      break;
      case 4 :
      Convole1DX_4( src,  r ,   a );
      break;
      case 5 :
      Convole1DX_5(src,  r ,   a );
      break;
      case 6 :
      Convole1DX_6(src,  r ,   a );
      break;
      case 7 :
      Convole1DX_7(src,  r,   a );
      break;
      case 8 :
      Convole1DX_8(src,  r ,   a );
      break;
      case 10 :
      Convole1DX_10(src,  r ,   a );
      break;
      case 11 :
      Convole1DX_11(src,  r ,   a );
      break;
      case 12 :
      Convole1DX_12(src,  r ,   a );
      break;
      case 13 :
      Convole1DX_13(src,  r ,   a );
      break;
      case 14 :
      Convole1DX_14(src,  r ,   a );
      break;
      case 15 :
      Convole1DX_15(src,  r ,   a );
      break;
      case 20 :
      Convole1DX_20(src,  r ,   a );
      break;
      default :
      Convole1DX_20(src,  r ,   a );

     }
  }
   CATCHALL{THROW_SAME;}
   ENDTRY

}

void Convole1DX(Image const& src, double r, int l , 
 Image const& mask, Image& a )
{ 
  int n = Demi_Fenetre(l,r);
  TRY{
    switch (n)
     {
      case 2 :
      Convole1DX_2( src,  r, mask,  a );
      break;
      case 3 :
      Convole1DX_3(src,  r, mask,  a );
      break;
      case 4 :
      Convole1DX_4(src,  r, mask,  a );
      break;
      case 5 :
      Convole1DX_5(src,  r, mask,  a );
      break;
      case 6 :
      Convole1DX_6(src,  r, mask,  a );
      break;
      case 7 :
      Convole1DX_7(src,  r, mask,  a );
      break;
      case 8 :
      Convole1DX_8(src,  r, mask,  a );
      break;
      case 10 :
      Convole1DX_10(src,  r, mask,  a );
      break;
      case 11 :
      Convole1DX_11(src,  r, mask,  a );
      break;
      case 12 :
      Convole1DX_12(src,  r, mask,  a );
      break;
      case 13 :
      Convole1DX_13(src,  r, mask,  a );
      break;
      case 14 :
      Convole1DX_14(src,  r, mask,  a );
      break;
      case 15 :
      Convole1DX_15(src,  r, mask,  a );
      break;
      case 20 :
      Convole1DX_20(src,  r, mask,  a );
      break;
      default :
      Convole1DX_20(src,  r,   mask,  a );

     }
    }
    CATCHALL{THROW_SAME;}
    ENDTRY

}

void Convole1DY (Image const& src, double r,  int l , Image& a )
{ 
  int n = Demi_Fenetre(l,r);

  TRY{
  switch (n)
     {
      case 2 :
      Convole1DY_2( src,  r ,    a );
      break;
      case 3 :
      Convole1DY_3( src,  r,    a );
      break;
      case 4 :
      Convole1DY_4( src,  r ,   a );
      break;
      case 5 :
      Convole1DY_5(src,  r ,   a );
      break;
      case 6 :
      Convole1DY_6(src,  r ,   a );
      break;
      case 7 :
      Convole1DY_7(src,  r,   a );
      break;
      case 8 :
      Convole1DY_8(src,  r ,   a );
      break;
      case 10 :
      Convole1DY_10(src,  r ,   a );
      break;
      case 11 :
      Convole1DY_11(src,  r ,   a );
      break;
      case 12 :
      Convole1DY_12(src,  r ,   a );
      break;
      case 13 :
      Convole1DY_13(src,  r ,   a );
      break;
      case 14 :
      Convole1DY_14(src,  r ,   a );
      break;
      case 15 :
      Convole1DY_15(src,  r ,   a );
      break;
      case 20 :
      Convole1DY_20(src,  r ,   a );
      break;
      default :
      Convole1DY_20(src,  r ,   a );

     }
  }
  CATCHALL{THROW_SAME;}
  ENDTRY

}

void Convole1DY(Image const& src, double r, int l , 
 Image const& mask, Image& a )
{
  int n = Demi_Fenetre(l,r);


  TRY{
 switch (n)
     {
      case 2 :
      Convole1DY_2( src,  r, mask,  a );
      break;
      case 3 :
      Convole1DY_3(src,  r, mask,  a );
      break;
      case 4 :
      Convole1DY_4(src,  r, mask,  a );
      break;
      case 5 :
      Convole1DY_5(src,  r, mask,  a );
      break;
      case 6 :
      Convole1DY_6(src,  r, mask,  a );
      break;
      case 7 :
      Convole1DY_7(src,  r, mask,  a );
      break;
      case 8 :
      Convole1DY_8(src,  r, mask,  a );
      break;
      case 10 :
      Convole1DY_10(src,  r, mask,  a );
      break;
      case 11 :
      Convole1DY_11(src,  r, mask,  a );
      break;
      case 12 :
      Convole1DY_12(src,  r, mask,  a );
      break;
      case 13 :
      Convole1DY_13(src,  r, mask,  a );
      break;
      case 14 :
      Convole1DY_14(src,  r, mask,  a );
      break;
      case 15 :
      Convole1DY_15(src,  r, mask,  a );
      break;
      case 20 :
      Convole1DY_20(src,  r, mask,  a );
      break;
      default :
      Convole1DY_20(src,  r,   mask,  a );

     }
 }
 CATCHALL{THROW_SAME;}
 ENDTRY

}

int
Demi_Fenetre_Optimale(double sigma,double precision)
{
  double demi = sigma*sqrt(-2. * log(precision)) ;
  int result = (int) ceil(demi) ;
  return(result) ;

}

//! returns the kernel actually used
Image ConvoleGauss1D1D( const Image &In, const double SigKernX, 
		       const double SigKernY, const double Precision,
		       Image &Out)
{
  Image imgcv1(In.Nx(), In.Ny());
  int lx = Demi_Fenetre(-1 ,SigKernX, Precision);
  Convole1DX(In, SigKernX, lx , imgcv1 ) ;
  double *kernx = Convole_MakeNoyau_1D(SigKernX, lx);

  int ly = Demi_Fenetre(-1 ,SigKernY, Precision);  
  Convole1DY(imgcv1, SigKernY, ly , Out);

  double *kerny = Convole_MakeNoyau_1D(SigKernY, ly);
  // in principle, the kernel normalization has been carried out.
  int nx = 2*lx+1;
  int ny = 2*ly+1;
  Image kern(nx,ny);
  for (int j=0;j<ny; ++j)
  for (int i=0; i<nx; ++i)
    kern(i,j) = kernx[i]*kerny[j];
  kern *= (1./kern.SumPixels());
  delete kernx; delete kerny;
  return kern;
}

Image ConvoleGauss1D1D( Image &InOut, const double SigKernX, 
		       const double SigKernY, const double Precision)
{
  return ConvoleGauss1D1D(InOut, SigKernX, SigKernY, Precision, InOut);
}

#ifdef FUTURE


/* Here are 1D convolution routines, somehow tested, in which the
   innermost loop runs over image pixels. This allows to optimize by
   hand as done in the following example:

static double scal_prod(double *x, double *y, const int size)
{
int nblock = size/10;
int remainder = size - nblock*10;
double sum = 0;
for (int i=0; i<nblock; ++i)
  {
  DO10( sum+= (*x)*(*y); ++x; ++y;)
  }
 for (int i=0; i<remainder; ++i) {sum+= (*x)*(*y); ++x; ++y;}
return sum;
}


These routine are already fast, and the scheme can easily be adapted to
2D convolution.

*/

bool ImConv_1DX(const Image &In, const double *K, const int S, Image &Out)
{
  int hs = S/2;
  if (2*hs+1 != S)
    {
      cout << " S should be odd ! " << endl;
      return false;
    }
  int nx = In.Nx();
  int ny = In.Ny();
  if (!(In.SameSize(Out)))
    {
      cout << " In aNd Out should have the same size !! " << endl;
      return false;
    }
  Out = 0;
  const int loops = nx-2*hs;

  for (int j=0; j<ny; ++j)
    {
      for (int k=-hs; k<=hs; ++k)
	{
	  Pixel *out = &Out(hs,j);
	  Pixel *in = &In(hs+k,j); // sign convention choosen here
	  double kernVal = K[k+hs];

	  // this loop to "unroll"
	  for (int kk=loops; kk; --kk)
	    {
	      *out += kernVal * (*in);
	      out++; in++;
	    }
	}
    }
  return true;
}


bool ImConv_1DY(const Image &In, const double *K, const int S, Image &Out)
{
  int hs = S/2;
  if (2*hs+1 != S)
    {
      cout << " S should be odd ! " << endl;
      return false;
    }
  int nx = In.Nx();
  int ny = In.Ny();
  if (!(In.SameSize(Out)))
    {
      cout << " In aNd Out should have the same size !! " << endl;
      return false;
    }
  Out = 0;
  const int loops = nx-2*hs;

  for (int j=hs; j<ny-hs; ++j)
    {
      for (int k=-hs; k<=hs; ++k)
	{
	  Pixel *out = &Out(0,j);
	  Pixel *in = &In(0,j+k); // sign convention choosen here
	  double kernVal = K[k+hs];

	  // this loop to "unroll"
	  for (int kk=loops; kk; --kk)
	    {
	      *out += kernVal * (*in);
	      out++; in++;
	    }
	}
    }
  return true;
}


bool EnlargeMask(Image &InOut, const int Border)
{

  int s = 2*Border+1;
  double *kernel1d = new double[s];
  for (int i=0; i<s; ++i) kernel1d[i] = 1;
  Image tmp(InOut.Nx(), InOut.Ny());;
  bool ok= (ImConv_1DX(InOut, kernel1d, s, tmp) && ImConv_1DY(tmp, kernel1d, s, InOut));
  delete [] kernel1d;
  return ok;
}


#endif /*FUTURE */

