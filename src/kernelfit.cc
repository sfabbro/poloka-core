#include <math.h>
#include <iostream>
#include <iomanip>

#include "kernelfit.h"
#include "basestar.h"
#include "vutils.h"

//#define DEBUG

#define DO2(A) A;A;
#define DO5(A) A;A;A;A;A;
#define DO10(A) DO2(DO5(A))

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


class Mat
{
  double *data;
  int d1,d2;
public :
  Mat(const int D1, const int D2) : d1(D1),  d2(D2)
  { data = new double[d1*d2]; memset(data,0,sizeof(double)*d1*d2); };
  ~Mat() { delete [] data;}
  double operator () (const int i, const int j) const { return data[i*d1+j];}
  double& operator () (const int i, const int j) { return data[i*d1+j];}
};

class Vect
{
  double *data;
  int d;
public :
  Vect(const int D): d(D)
  { data = new double[d]; memset(data,0,sizeof(double)*d); };
  ~Vect() { delete [] data;}
  double operator () (const int i) const { return data[i];}
  double& operator () (const int i) { return data[i];}
};


int KernelFit::FitDifferentialBackground(const double NSig)
{  
   
  if (optParams.SepBackVar.Degree == -1) return 0;

  unsigned int nterms = optParams.SepBackVar.Nterms();

  diffbackground.resize(nterms);
  WorstDiffBkgrdSubtracted = new Image(WorstImage->Nx(),WorstImage->Ny()); 
  

  Mat A(nterms,nterms);
  Vect B(nterms);
  Vect monom(nterms);

  Pixel bestMean,bestSig,worstMean,worstSig;
  BestImage->SkyLevel(DataFrame,&bestMean, &bestSig);
  WorstImage->SkyLevel(DataFrame,&worstMean, &worstSig);
  
  double cut1  = NSig*bestSig;
  double cut2  = NSig*worstSig;

  int ibeg,iend,jbeg,jend;
  ibeg = int(DataFrame.xMin);
  iend = int(DataFrame.xMax);
  jbeg = int(DataFrame.yMin);
  jend = int(DataFrame.yMax);

  int jump = 10;
  for (int j=jbeg ; j< jend ; j+= jump)
  for (int i=ibeg ; i< iend ; i+= jump)
    {
      Pixel p1 = (*WorstImage)(i,j);
      if (fabs(p1-bestMean) > cut1) continue;
      Pixel p2 = (*BestImage)(i,j);
      if (fabs(p2-worstMean) > cut2) continue;
      for (int q1=0; q1<nterms; ++q1)
	{
	  monom(q1) = optParams.SepBackVar.Value(double(i), double(j),q1);

	  for (int q2 = q1; q2<nterms; ++q2) A(q1,q2) += monom(q1)*monom(q2);
	  B(q1) += monom(q1)*p1;
	}
    }
  /* symetrize */
  for (int q1=0; q1<nterms; ++q1) 
    for (int q2 = q1+1; q2<nterms; ++q2) A(q2,q1) = A(q1,q2); 
  if (!MatSolve(&A(0,0),nterms,&B(0)))
    {
      cerr << " could not compute differential background !!!!" << endl;
      return 0;
    }
  
  cout << setprecision(10);
  cout << " separately fitted differential background " << endl;
  cout << " ----------------------------------------- " << endl;
  for (int q1=0; q1< nterms; ++q1) cout << B(q1) << " " ;
  cout << endl;

  for (int q1=0; q1<nterms; ++q1)
      diffbackground[q1] = B(q1);

  for (int j=jbeg ; j< jend ; ++j)
  for (int i=ibeg ; i< iend ; ++i)
    {
      (*WorstDiffBkgrdSubtracted)(i,j) = (*WorstImage)(i,j) - BackValue(i,j);
    }
  return 1;
}



#define OPTIMIZED /* means pushing pointers by hand ... */

// #define DO3(I) I;I;I;
// #define DO9(I) DO3(DO3(I))
// #define DO19(I) DO9(I);DO9(I);I

#include <time.h>
//static double sqr(double x) { return x*x;}

//! convolves the best image with the current kernel (Usually set by DoTheFit).
/*! The UpdateKern parameter is the pixel range over which the kernel will not
be updated. The concolved image is to be found in ConvolvedBest.*/

void KernelFit::BestImageConvolve(int UpdateKern)
{
clock_t tstart = clock();
const Image &Source = *BestImage;
// copy the input image so that side bands are filled with input values.
ConvolvedBest = new Image(Source); 
if (solution.size()==0)
   {
     cerr << "BestImageConvolve : the kernel fit was not done yet... or was impossible .. no convolution" << endl;
     return;
   }
int ksx = optParams.HKernelSize;
int ksy = optParams.HKernelSize;
int startx = ksx;
int starty = ksy;
int endx = Source.Nx() - ksx;
int endy = Source.Ny() - ksy;
Kernel kern(ksx,ksy);
int npix = kern.Nx();
int nxregions = (UpdateKern) ? (Source.Nx()/UpdateKern)+1 : 1;
int nyregions = (UpdateKern) ? (Source.Ny()/UpdateKern)+1 : 1;
cout << " Convolving best image" << endl;
// some printout:
KernCompute(kern, Source.Nx()/2, Source.Ny()/2);
 {cout <<  " kernel caracteristics at i j " <<  Source.Nx()/2 << ' ' << Source.Ny()/2; kern.dump_info();}
 for (int i = 1 ; i <= 3 ; i++)
   for (int j = 1 ; j <= 3 ; j++)
     {
       KernCompute(kern, i*Source.Nx()/4, j*Source.Ny()/4);
       cout <<  " kernel caracteristics at i j " <<  i*Source.Nx()/4 
	    << ' ' << j*Source.Ny()/4; kern.dump_info();
     }

 // on the road again
for (int iry = 0; iry < nyregions; ++iry)
for (int irx = 0; irx < nxregions; ++irx)
  {
  int sj = max(iry*UpdateKern,starty);
  int ej = min(sj+UpdateKern,endy);
  int si = max(irx*UpdateKern, startx);
  int ei = min(si+UpdateKern, endx);
  double xc = (si+ei)/2;
  double yc = (sj+ej)/2;
  KernCompute(kern, xc, yc);
  double kernSum = kern.sum();

  if (((irx == 0) || (irx == nxregions-1)) &&  ((iry ==0) || (iry == nyregions - 1)))
    {cout <<   " kernel caracteristics at i j " << xc << ' ' << yc; kern.dump_info();}
  for (int j = sj; j < ej; ++j)
  for (int i = si; i < ei; ++i)
    {
      double sum = 0;
#ifndef OPTIMIZED
      for (int jk = -ksy; jk <=ksy; ++jk)
      for (int ik = -ksx; ik <=ksx; ++ik)
         sum += kern(ik,jk)*Source(i-ik,j-jk);
#else
      DPixel *pk = kern.begin();
      for (int jk = -ksy; jk <= ksy; ++jk)
	{
	  Pixel *ps = &Source(i+ksx, j-jk);
 	  for (int k=npix; k; --k) {sum += (*pk) * (*ps); ++pk ; --ps;}
	}
#endif
      sum -= kernSum*BestImageBack;
  //  if (kern.begin() + kern.Nx()*kern.Ny() - pk ) cout << " BestImageImageConvol catastrophe...." << endl;
      (*ConvolvedBest)(i,j) = sum;
    }
  }
/* account for differential background. We do it for the whole image, including
   side bands where the convolution did not go but where anyway initialized
   to the input value. TODO : correct side bands for photometric ratio */

 int sx = ConvolvedBest->Nx();
 int sy = ConvolvedBest->Ny();
 for (int j=0; j < sy; ++j) for (int i=0; i < sx ; ++i)
  (*ConvolvedBest)(i,j) += BackValue(i,j) + WorstImageBack;

clock_t tend = clock();
cout << " CPU for convolution " << float(tend-tstart)/float(CLOCKS_PER_SEC) << endl;
}

Image *KernelFit::VarianceConvolve(const Image &Source, int UpdateKern)
{
clock_t tstart = clock();
Image *Result = new Image(Source);
if (solution.size()==0)
   {
     cerr << "Variance Convolve : the kernel fit was not done yet... or was impossible .. no convolution" << endl;
     return NULL;
   }
int ksx = optParams.HKernelSize;
int ksy = optParams.HKernelSize;
int startx = ksx;
int starty = ksy;
int endx = Source.Nx() - ksx;
int endy = Source.Ny() - ksy;
Kernel kern(ksx,ksy);
int npix = kern.Nx();
int nxregions = (UpdateKern) ? (Source.Nx()/UpdateKern)+1 : 1;
int nyregions = (UpdateKern) ? (Source.Ny()/UpdateKern)+1 : 1;
cout << " Convolving variance" << endl;
// some printout:
KernCompute(kern, Source.Nx()/2, Source.Ny()/2);
 {cout <<  " kernel caracteristics at i j " <<  Source.Nx()/2 << ' ' << Source.Ny()/2; kern.dump_info();}
 for (int i = 1 ; i <= 3 ; i++)
   for (int j = 1 ; j <= 3 ; j++)
     {
       KernCompute(kern, i*Source.Nx()/4, j*Source.Ny()/4);
       kern *= kern;
       cout <<  " kernel^2 caracteristics at i j " <<  i*Source.Nx()/4 
	    << ' ' << j*Source.Ny()/4; kern.dump_info();
     }

 // on the road again
for (int iry = 0; iry < nyregions; ++iry)
for (int irx = 0; irx < nxregions; ++irx)
  {
  int sj = max(iry*UpdateKern,starty);
  int ej = min(sj+UpdateKern,endy);
  int si = max(irx*UpdateKern, startx);
  int ei = min(si+UpdateKern, endx);
  double xc = (si+ei)/2;
  double yc = (sj+ej)/2;
  KernCompute(kern, xc, yc);
  
  // This is a variance :
  kern *= kern;

  for (int j = sj; j < ej; ++j)
  for (int i = si; i < ei; ++i)
    {
      double sum = 0;
#ifndef OPTIMIZED
      for (int jk = -ksy; jk <=ksy; ++jk)
      for (int ik = -ksx; ik <=ksx; ++ik)
         sum += kern(ik,jk)*Source(i-ik,j-jk);
#else
      DPixel *pk = kern.begin();
      for (int jk = -ksy; jk <= ksy; ++jk)
	{
	  Pixel *ps = &Source(i+ksx, j-jk);
 	  for (int k=npix; k; --k) {sum += (*pk) * (*ps); ++pk ; --ps;}
	}
#endif
      (*Result)(i,j) = sum;
    }
  }
clock_t tend = clock();
cout << " CPU for convolution " << float(tend-tstart)/float(CLOCKS_PER_SEC) << endl;
 return Result;
}


static double my_pow(const double X, const int k)
{
  if (k==0) return 1;
  double res = X;
  for (int i=1; i<k; ++i) res *= X;
  return res;
}


static void PolGaussKern(Kernel &Kern, double Sig, int DegX, int DegY)
{
  /* fills the kernel array with values : exp (-(x**2+y**2)/(2.sigma**2))*(x**degx)*(y**degy) */
int s = Kern.HSizeX();
int totSize = 2*s+1;
double *vx = new double[totSize];
double *vy = new double[totSize];
double alpha = 1./(2.*Sig*Sig);
for (int i=0; i<totSize; ++i)
  {
  double x = double(i-s);
  double valg = exp(-x*x*alpha);
  vx[i] = valg*my_pow(x,DegX);
  vy[i] = valg*my_pow(x,DegY);
  }

DPixel *p = Kern.begin();
for (int j=0;j<totSize;++j)
  for (int i=0; i<totSize; ++i)
    {*p = vx[i]*vy[j]; ++p;}
delete [] vx;
delete [] vy;
}
 

static void SetDelta(Kernel &Kern)
{
Kern.Zero();
Kern(0,0) = 1.0;
}



void alloc_m_and_b(vector<double>  &m, vector<double> &b, int mSize)
{
  m.resize(mSize*mSize);
  b.resize(mSize*mSize);
}

void KernelFit::KernelsFill()
{

if (Kernels.size() != 0)
  {
    cout << "entering KernelsFill with a non empty Kernel array ... " << endl;
    Kernels.clear(); 
  }
/* all kernels should have the same size */
Kernels.push_back(Kernel(optParams.HKernelSize, optParams.HKernelSize));
SetDelta(Kernels[0]);
 int oldprec = cout.precision();
 cout << setprecision(10);
cout << "  kernels parameters (sigma, degree) : ";
for (int iwidth = 0; iwidth < optParams.NGauss ; ++iwidth) 
       cout << '(' <<optParams.Sigmas[iwidth] << ',' << optParams.Degrees[iwidth] << ')';
cout << endl;  

for (int iwidth = 0; iwidth < optParams.NGauss ; ++iwidth)
  {
  int maxdeg = optParams.Degrees[iwidth];
  for (int degx = 0; degx <= maxdeg; ++degx)
    for (int degy =0; degy<= maxdeg-degx; ++degy)
      {
	Kernel kern(optParams.HKernelSize, optParams.HKernelSize);
	//        cout << " degx degy " << degx << ' ' << degy << endl;
	PolGaussKern(kern,optParams.Sigmas[iwidth], degx, degy);
        Kernels.push_back(kern);
      }
  }
mSize = Kernels.size() * optParams.KernVar.Nterms() + optParams.BackVar.Nterms();
 cout << " KernelsFill : Number of basic kernels : " << Kernels.size()
      << " number of fitted coefficients "  << mSize <<endl;

 cout << setprecision(oldprec);
}



static double image_scal_prod(const DImage &Vignette, const Image* I, const double IBack, int xs, int ys)
{
double sum = 0;
int xsize = Vignette.Nx();
int ysize = Vignette.Ny();
for (int j=0; j< ysize; ++j)
  {
  for (int i=0; i< xsize; ++i)
    {
     
    sum += (double((*I)(i+xs,j+ys)) - IBack)* Vignette(i,j);
    } 
  }
return sum;
}


static double image_sum(const Image* I, const double IBack, int xs, int ys, int xe, int ye)
{
double sum = 0;
for (int j=ys; j<=ye ; ++j)
  {
  for (int i=xs; i<=xe; ++i)
    {
    sum += double((*I)(i,j)) - IBack;
    } 
  }
return sum;
}

void KernelFit::AllocateConvolvedStamps() /* allocate a C-vector of DImages to store a_stamp x Kernels[i] */
{
if (convolutions) return;
int nkern = Kernels.size();
convolutions = new DImage[nkern];
int convolvedSize = optParams.ConvolvedSize();
for (int i=0; i<nkern; ++i) 
  convolutions[i].Allocate(convolvedSize,convolvedSize);
int nBackStamps = optParams.BackVar.Nterms();
backStamps = new DImage[nBackStamps];
for (int i=0; i<nBackStamps; ++i) 
  backStamps[i].Allocate(convolvedSize,convolvedSize);
}

void KernelFit::DeallocateConvolvedStamps()
{
if (convolutions) 
  {
    delete [] convolutions;
    convolutions = NULL;
    delete [] backStamps;
    backStamps = NULL;
  }
}


void KernelFit::OneStampMAndB(const Stamp &AStamp, vector<double> &StampM, vector<double> &StampB)
{
  int oldprec = cout.precision();
 cout << setprecision(10);
 ios::fmtflags  old_flags = cout.flags(); 
 cout << resetiosflags(ios::fixed) ;
 cout << setiosflags(ios::scientific) ;

 int nkern = Kernels.size();
/* assume that all kernels have the same size: */
 //int ksize = optParams.HKernelSize;
/* assume that all stamps have the same size: */
 //int stampsize = AStamp.Nx();
 int convolvedSize = optParams.ConvolvedSize(); 

 /* 'convolutions' and backStamps are parts of KernelFit that are allocated 
    by AllocateConvolvedStamps and deallocated by Dealloc... done this 
    way so that it is allocated once in the loop on stamps :*/

 AllocateConvolvedStamps(); /* does nothing if already done */

 int hConvolvedSize = convolvedSize/2;
 int convolvedPix = convolvedSize*convolvedSize;
 double *spatialCoeff = new double [optParams.KernVar.Nterms()];
 // double *backCoeff = new double [optParams.BackVar.Nterms()];

 /* values of monomials such as (xc**n)*(yc**m), for this stamp  : */
 // for the smooth kernel variation
 for (int i = 0; i  < optParams.KernVar.Nterms(); ++i) 
     spatialCoeff[i] = optParams.KernVar.Value(AStamp.xc, AStamp.yc, i);
 //for the background, actually compute monomials for all stamp pixels:
 for (int ib =0; ib < optParams.BackVar.Nterms(); ++ib) 
   {
     DImage &backStamp = backStamps[ib];
     for (int j= 0; j < convolvedSize ; ++j)
     for (int i= 0; i < convolvedSize ; ++i)
       {
	 int xi = i+ AStamp.xc - hConvolvedSize;
	 int yj = j+ AStamp.yc - hConvolvedSize;
	 backStamp(i, j) = optParams.BackVar.Value(xi, yj, ib);
       }
   }

 StampM.resize(mSize*mSize);
 StampB.resize(mSize);

    /* contributions to the matrix m */
    /* background-background terms */
 for (int ib = 0; ib < optParams.BackVar.Nterms(); ++ib)
   {
   for (int jb =0; jb<=ib; ++jb)   
     {
     StampM[BackIndex(ib)*mSize+BackIndex(jb)] += 
       scal_prod(backStamps[ib].begin(), backStamps[jb].begin(), convolvedPix);
     }
   }
 
 for (int ik=0; ik < nkern; ++ik)
   {
    /* kernel - kernel terms */

   Convolve(convolutions[ik], AStamp.source, Kernels[ik]);

   for (int jk=0; jk <=ik; ++jk)
     {
     double integral = scal_prod(convolutions[ik].begin(), 
				 convolutions[jk].begin(),
				 convolvedPix);
     for (int is =0; is < optParams.KernVar.Nterms(); ++is)
       {
       int im = KernIndex(ik,is);
       for (int js =0; js < optParams.KernVar.Nterms(); ++js)
	 {
	 int jm = KernIndex(jk,js);
         StampM[im*mSize+jm] += integral*spatialCoeff[is]*spatialCoeff[js];  // alard eq (3)
         }
       }
     }

      /* kernel-background terms */
   for (int is = 0; is < optParams.KernVar.Nterms(); ++is)
     {
     int im = KernIndex(ik,is);
     for (int jb=0; jb < optParams.BackVar.Nterms(); ++jb) 
       {
       int jm = BackIndex(jb);
       /* fill only m part for j<=i 
	 StampM(i,j) is in principle [i*mSize+j], swap them here...*/
       StampM[jm*mSize+im] += spatialCoeff[is]*
	 scal_prod(convolutions[ik].begin(), backStamps[jb].begin(), 
		   convolvedPix);
       }
     }

      /* compute contributions to b */

      /* kernel terms */
   const Image *worst;
   if (optParams.SepBackVar.Nterms() != 0)
     {
       worst = WorstDiffBkgrdSubtracted;
     }
   else
     {
      worst = WorstImage;
     }
   
   double bintegral = image_scal_prod(convolutions[ik], worst, WorstImageBack, AStamp.xc - hConvolvedSize, AStamp.yc - hConvolvedSize);
   for (int is =0; is < optParams.KernVar.Nterms(); ++is)
     {
     int im = KernIndex(ik,is);
     StampB[im] += bintegral * spatialCoeff[is]; // alard eq (4)
     }
   } /* end of for (ik */

      /* background terms of b */
 for (int ib=0; ib < optParams.BackVar.Nterms(); ++ib)
   {
   StampB[BackIndex(ib)] += 
     image_scal_prod(backStamps[ib], WorstImage, WorstImageBack, 
		     AStamp.xc - hConvolvedSize, AStamp.yc - hConvolvedSize);
   }
 delete [] spatialCoeff;
 // delete [] backCoeff;
 /* return a symtetrized matrix ... safer than too assume anything 
    about calling routine*/
for (int i=0; i<mSize; ++i) for (int j=i+1; j<mSize; ++j) StampM[i*mSize+j] = StampM[j*mSize+i];

 cout << setprecision(oldprec);
 cout.flags(old_flags);
}


void KernelFit::SubtractStampFromMAndB(Stamp& AStamp)
{
 vector<double> mStamp;
 vector<double> bStamp;
 alloc_m_and_b(mStamp, bStamp, mSize);
 OneStampMAndB(AStamp, mStamp, bStamp);
 for (int i = mSize*mSize-1; i>=0; --i) m[i] -= mStamp[i];
 for (int i = mSize-1; i>=0; --i) b[i] -= bStamp[i];
}



void KernelFit::ComputeMAndB()
{

/* allocate the convolved stamps */


alloc_m_and_b(m,b,mSize);

 vector<double> mStamp;
 vector<double> bStamp;
alloc_m_and_b(mStamp, bStamp, mSize);


for (StampIterator si = BestImageStamps->begin(); si != BestImageStamps->end(); ++si)
  {
    /* since the allocation of convolved stamps once for all assumes that all stamps have the same size, 
       should check it */
    OneStampMAndB(*si, mStamp, bStamp);
    for (int i = mSize*mSize-1; i>=0; --i) m[i] += mStamp[i];
    for (int i = mSize-1; i>=0; --i) b[i] += bStamp[i];

  } /* end of loop on stamps */

#ifdef DEBUG
 int oldprec = cout.precision();
 cout << setprecision(10);
 ios::fmtflags  old_flags = cout.flags(); 
 cout << resetiosflags(ios::fixed) ;
 cout << setiosflags(ios::scientific) ;
 cout << " matrix " << endl;
 for (int i=0; i<mSize; ++i) { for (int j=0; j<=i; ++j) cout << m[i*mSize+j] << ' '; cout << endl;}
 cout << " b " << endl;
 for (int i=0; i<mSize; ++i) cout << b[i] << ' ' ; cout <<endl;
 cout << setprecision(oldprec);
 cout.flags(old_flags);
#endif


} 


static void KernLinComb(Kernel &Result, const vector<Kernel> &VK, double *Coeffs)
{
int size = Result.Nx() * Result.Ny();
int nkern = int(VK.size());
Result.Zero();
for (int ik=0; ik < nkern; ++ik)
  {
  DPixel *k = VK[ik].begin();
  DPixel *r = Result.begin();
  double coeff = Coeffs[ik];
  for (int i=0; i< size; ++i) {*r += *k * coeff; ++r ; ++k;}
  }
}


#include "toadscards.h"
#include "datacards.h"
#include "fileutils.h" // for FileExists

OptParams::OptParams()
{
NGauss = 3;
Sigmas.resize(NGauss);
Degrees.resize(NGauss);
Sigmas[0] = 0.7; Sigmas[1] = 1.5; Sigmas[2] = 2.;
Degrees[0] = 6; Degrees[1] = 4; Degrees[2] = 2;
//Degrees[0] = 1; Degrees[1] = 4; Degrees[2] = 2;
HStampSize = 15;
HKernelSize = 9;
MaxStamps = 150;
KernVar.SetDegree(2); // degree of spatial variations of the kernel 
BackVar.SetDegree(2); //degree of spatial variations of the background
SepBackVar.SetDegree(-1); //degree of spatial variations of the background if you want to fit it separately
 NSig = 4;
 string dataCardsName = DefaultDatacards();
 if (FileExists(dataCardsName))
   {
     DataCards cards(dataCardsName);
     if (cards.HasKey("KFIT_SIG_GAUSS"))
       {
	 int nGauss = cards.NbParam("KFIT_SIG_GAUSS");
	 int nGaussDeg = cards.NbParam("KFIT_DEG_GAUSS");
	 if (nGaussDeg != nGauss)
	   {
	     cerr << dataCardsName 
		  << " : don't find the same number of items in KFIT_SIG_GAUSS and KFIT_DEG_GAUSS" 
		  << endl;
	     cerr << " cutting the longest one  to the size of the shortest one "<< endl;
	     nGauss = min(nGauss,nGaussDeg);
	   }
	 if (nGauss > 0) 
	   {
	     NGauss = nGauss;
	     Sigmas.resize(NGauss);
	     Degrees.resize(NGauss);
	     for (int i=0; i<NGauss; ++i)
	       {
		 Sigmas[i] = cards.DParam("KFIT_SIG_GAUSS", i);
		 Degrees[i] = cards.IParam("KFIT_DEG_GAUSS",i);
	       }	   
	   }
       } // HasKey("KFIT_SIG_GAUSS")
     if (cards.HasKey("KFIT_NSIG")) NSig = cards.DParam("KFIT_NSIG");
     if (cards.HasKey("KFIT_KERNVAR_DEG")) 
       {
	 int deg = cards.IParam("KFIT_KERNVAR_DEG");
	 KernVar.SetDegree(deg); // degree of spatial variations of the kernel 
       }
     if (cards.HasKey("KFIT_BACKVAR_DEG")) 
       {
	 int deg = cards.IParam("KFIT_BACKVAR_DEG");
	 BackVar.SetDegree(deg); // degree of spatial variations of the background
       }
     if (cards.HasKey("KFIT_SEPBACKVAR_DEG")) 
       {
	 int deg = cards.IParam("KFIT_SEPBACKVAR_DEG");
	 SepBackVar.SetDegree(deg); // degree of spatial variations of the background
       }
   } // if (has datacards)
}
  

void OptParams::OptimizeSizes(double BestSeeing, double WorstSeeing)
{
 int oldprec = cout.precision();
 cout << setprecision(10);
cout << " Choosing sizes with bestseeing = " <<  BestSeeing << " WorstSeeing = " << WorstSeeing << endl;
if ( BestSeeing > WorstSeeing) swap(BestSeeing, WorstSeeing);
double kernSig = max(sqrt(WorstSeeing*WorstSeeing - BestSeeing*BestSeeing),0.4);
 cout << " Expected kernel sigma  : " << kernSig << endl ; 
HKernelSize = max(int( ceil (NSig * kernSig)),4);
 for (int i=0; i<NGauss; ++i) Sigmas[i] *= kernSig;
/* convolvedStampSize is the size of the stamps that enter the chi2. 2 is added
to get some area to estimate the differential background. */
int convolvedStampSize = int(ceil(NSig * WorstSeeing)) + 2;
HStampSize = HKernelSize + convolvedStampSize;
 cout << setprecision(oldprec);
}

void OptParams::OptimizeSpatialVariations(const int NumberOfStars)
{
  int degree = KernVar.Degree;
  if (NumberOfStars < 75) degree--;
  if (NumberOfStars < 50) degree--;
  // only lower greees, not raise it from datacards value
  if (KernVar.Degree > degree) 
    {
      cout << " Lowering spatial kernel (deg = " << degree 
	   << ") variations with NumberOfStars = " << NumberOfStars << endl;
      KernVar.SetDegree(degree);
    }
  //  if (BackVar.Degree > degree) BackVar.SetDegree(degree); not so useful

  if (NumberOfStars < 25)
    {
      cout << " Choosing spatial variations with NumberOfStars = " 
	   << NumberOfStars << endl;
      cout << " lowered degree of polynomials associated to gaussians to ";
      for (int i=0; i<NGauss; ++i) {Degrees[i] /= 2; cout << Degrees[i] << ' ';}
      cout << endl;
    }

  if (NumberOfStars < 8) 
    {
      cout << " lowered degree of polynomials associated to gaussians to ";
      for (int i=0; i<NGauss; ++i) {Degrees[i] /= 2; cout << Degrees[i] << ' ';}
      cout << endl;
    }
}

void XYPower::SetDegree(const int DegreeValue)
{
  Degree = DegreeValue;
  int N_terms = (Degree+1)*(Degree+2)/2;
  Xdeg.resize(N_terms);
  Ydeg.resize(N_terms);
  int q = 0;
  for (int xdeg=0; xdeg <=Degree; ++xdeg)
    for (int ydeg = 0; ydeg <= Degree - xdeg; ++ydeg)
      {
	Xdeg[q] = xdeg;
	Ydeg[q] = ydeg;
	++q;
      }
}



double XYPower::Value(const double X, const double Y, const int q) const
{
  if (q>=Nterms()) {cerr << "  XYPower::Value ..."  << endl; abort();}
// my_pow(double,int) is about 5 times faster than pow(double,doublke)
  return my_pow(X/100.,Xdeg[q])*my_pow(Y/100.,Ydeg[q]); 
}


void KernelFit::KernCompute(Kernel &Kern, const double X, const double Y) const
{
  int nKern = Kernels.size();
  double *coeff = new double[nKern];
  for (int i=0; i<nKern; ++i) coeff[i] = 0;
  for (int is =0; is < optParams.KernVar.Nterms(); ++is)
    {
      double spatialCoeff = optParams.KernVar.Value(X,Y,is);
      for (int ik=0; ik< nKern; ++ik)
    {  
      coeff[ik] +=  solution[KernIndex(ik, is)]*spatialCoeff;
    }
    }
  KernLinComb(Kern, Kernels, coeff);
  delete [] coeff;
}

static void mean_median_sigma(double *values, const int nval, double &mean,  double &median, double &sigma)
{
  mean =0;
  sigma =0;
  median = DArrayMedian(values, nval);
  for (int i=0; i<nval-1 ; ++i)
    {
      mean += values[i];
      sigma += values[i]*values[i];
    }
  mean /= double(nval);
  sigma = sigma/double(nval) - mean*mean;
  if (sigma>0)  sigma = sqrt(sigma); else sigma = 0;
}

double KernelFit::StampChi2(Stamp &stamp, double VSky, double InvGain)
{
Kernel kern(optParams.HKernelSize, optParams.HKernelSize);
int  xc = stamp.xc;
int  yc = stamp.yc;
KernCompute(kern, xc, yc);
int convolvedSize = optParams.ConvolvedSize(); /* should be odd */
int hConvolvedSize = convolvedSize/2;
DImage R_conv_K(convolvedSize,convolvedSize);
Convolve(R_conv_K, stamp.source, kern);
double chi2 = 0;
int xs = xc - hConvolvedSize; 
int ys = yc - hConvolvedSize;
for (int j=0; j< convolvedSize; ++j)
for (int i=0; i< convolvedSize; ++i)
  {
  double w_value = (*WorstImage)(i+xs, j + ys);
  double res = R_conv_K(i,j) - w_value + WorstImageBack + BackValue(i + xs,j + ys); 
  // from the computation of chi2 , why chi2<0 are due to bad columns? (see above comment)
  // it looks like it is more due to a negative fluctuation.
  // furthermore, although we weight for chi2, we don't seem to weight 
  // for the matrix-vector filling.

  chi2 += res*res/(InvGain*w_value+VSky);
  }
stamp.chi2 = chi2;
return chi2;
}


/* removal of stamps which are found different after kernel fit, convolution and subtraction */
void KernelFit::FilterStamps()
 {
 int dropped;
 int nstamps =  BestImageStamps->size();
 double *chi2s = new double[nstamps];
 int oldprec = cout.precision();
 cout << setprecision(10);
 cout << " sigma(sky) of worst image " << sqrt(SkyVarianceWorstImage) << endl;
 cout << " BackValue(0,0)" << BackValue(0,0) << endl;
 do 
   {
     int i=0;
     for (StampIterator si = BestImageStamps->begin(); si != BestImageStamps->end(); ++si, ++i)
       {
	 Stamp &stamp = *si;
	 chi2s[i] =  StampChi2(stamp, SkyVarianceWorstImage, 1./WorstImageGain); /* also assigns stamp.chi2 */
	 // cout << " xc yc chi2 " << stamp.xc << ' ' << stamp.yc << ' ' << chi2s[i] << endl;
	 // if (stamp.star) stamp.star->dump();
       }  
     double mean,sigma,median;
     mean_median_sigma(chi2s,nstamps,mean,median,sigma);
     cout << " mean median sigma " << mean << ' ' << median << ' ' << sigma << endl;
     cout << " mean/dof " << mean/(optParams.ConvolvedSize()*optParams.ConvolvedSize()) << endl;
     dropped = 0;
     
     // this trick is used to fit the kernel using a predefined catalog
     // of objet to compute the kernel
     if (getenv("NOFILTERING")) break;
     for (StampIterator si = BestImageStamps->begin(); si != BestImageStamps->end(); )
       {
	 Stamp &stamp = *si; 
	 /* cut too large chi2's and negative ones (due to dead columns)*/
	 if (stamp.chi2 > median + 4*sigma || stamp.chi2 < 0) 
	   {
	     cout << " delete : xc yc chi2 " << stamp.xc << ' ' << stamp.yc << ' ' << stamp.chi2 << endl;
	     //if (stamp.star) stamp.star->dump();
	     dropped++;
	     SubtractStampFromMAndB(stamp);
	     si = BestImageStamps->erase(si);
	   }
	 else ++si;
       }
     if (dropped) /* we have to recompute the kernel (just solve once again the normal equation) */
       {
	 cout << " dropped " << dropped << " stamps : refitting " << endl;
	 Solve();
       }
     nstamps -= dropped;
   } 
 while (dropped != 0);
 
 cout  << " finished the fit  with " << BestImageStamps->size() << " stamps " << endl;
 if (getenv("DUMP_FIT_LIST"))
   {
     int i=0;
     for (StampIterator si = BestImageStamps->begin(); si != BestImageStamps->end(); ++si, ++i)
       {
	 cout << "Chi2 " << (*si).chi2 << " " << *(*si).star;
       }
   }
 double mean,median,sigma;
 mean_median_sigma(chi2s, nstamps, mean, median, sigma);
 chi2 = mean/(optParams.ConvolvedSize()*optParams.ConvolvedSize());
 cout << " final mean/dof " << chi2 << endl;
 cout << setprecision(oldprec);
  
delete [] chi2s;
}


#define CONST_INTEGRAL

int KernelFit::Solve()
{
 int oldprec = cout.precision();
 cout << setprecision(10);
 ios::fmtflags  old_flags = cout.flags(); 
 cout << resetiosflags(ios::fixed) ;
 cout << setiosflags(ios::scientific) ;
if (solution.size()==0) solution.resize(mSize);

#ifdef CONST_INTEGRAL /* i.e. the integral of the kernel (photometric ratio) is constant over the image */
// use Lagrange multipliers technique
 cerr <<"integral of kernel is assumed to be constant." << endl ;
int nKern = Kernels.size();
int nc = optParams.KernVar.Nterms() -1; // number of constraints
int totSize = mSize + nc;

//double *mprime = new double [totSize*totSize];
//memset(mprime,0,sizeof(double)*totSize*totSize);
//double *bprime = new double [totSize];
//memset(bprime,0,sizeof(double)*totSize);
//memcpy(bprime,b,sizeof(double)*mSize);
 
 vector<double> mprime;
 vector<double> bprime;
 mprime.resize(totSize*totSize);
 bprime = b;
 bprime.resize(totSize);
 

for (int i=0; i<mSize; ++i)
  for (int j=0; j<mSize; ++j)
    {
    mprime[i*totSize+j] = m[i*mSize+j];
    }
for (int ik=0; ik < nKern; ++ik)
  {
  double kern_int = Kernels[ik].sum();
  for (int ic =1; ic < optParams.KernVar.Nterms(); ++ic)
    {
    int ip = KernIndex(ik,ic);
    int jp = mSize+ic-1;
    mprime[ip*totSize+jp] = mprime[jp*totSize+ip] = kern_int;
    }
  }
// have to use a lin eq. solver that accomodates non posdef matrices : mprime is NOT posdef.
int inversion=MatSolve(&mprime[0],totSize,&bprime[0]);
cout << " inversion " << inversion << endl;
if (inversion)
  {
    solution = bprime;
  }

#else /* NO CONSTRAINT on the variations of kernel integral */
/* operate on a copy, to preserve m */

 vector<double> mprime;
 vector<double> bprime;
 mprime = m;
 bprime = b;
 

int inversion = MatSolve(mprime, mSize, bprime);
 cout << " inversion " << inversion << endl;
if (inversion)
  {
    solution = bprime;
  }
#endif

Kernel kernel_at_center( optParams.HKernelSize, optParams.HKernelSize);
KernCompute(kernel_at_center, BestImage->Nx()/2, BestImage->Ny()/2);

/* cout << "solution" << endl;
 for (int i=0; i<mSize; ++i) { 
   cout << solution[i] << ' ' ; 
   if ( (i%10) == 9) cout <<endl;}
   cout <<endl; */

 KernAtCenterSum = kernel_at_center.sum();
 cout << " Kernel integral ( = photometric ratio) " << KernAtCenterSum  << endl;
 if (optParams.BackVar.Nterms()>0) 
   {
     cout << setprecision(10);
     cout << " differential background " << endl;
     for (unsigned int ib=0; ib< optParams.BackVar.Nterms(); ++ib) cout << solution[BackIndex(ib)] << " " ;
     cout << endl;
   }

 double vx,vy,vxy; kernel_at_center.moments(vx,vy,vxy);
 cout << " sigmas " << sqrt(vx) << ' ' << sqrt(vy) << " rho " << vxy/sqrt(vx*vy) << endl;
 // cout << " Kernel " << endl; kernel_at_center.dump();
 cout << setprecision(oldprec);
 cout.flags(old_flags);
 return inversion;
}

double KernelFit::BackValue(const double&x, const double &y) const
{
  double val=0.0;

  if (optParams.BackVar.Nterms()==0)
    {
      val= diffbackground[0]; /* the constant in the polynomial .. */
      for (unsigned int ib=1; ib < optParams.SepBackVar.Nterms(); ++ib)
	{
	  val += diffbackground[ib]*optParams.SepBackVar.Value(x,y,ib);
	}
    }
  else
    {
      val= solution[BackIndex(0)]; /* the constant in the polynomial .. */
      for (unsigned int ib=1; ib < optParams.BackVar.Nterms(); ++ib)
	{
	  val += solution[BackIndex(ib)]*optParams.BackVar.Value(x,y,ib);
	}
    }
  return val;
}


//! does the kernel fit (DoTheFit) and the convolution (BestImageConvolve).
int KernelFit::DoIt(const BaseStarList &List, double &BestSeeing, double &WorstSeeing)
{
DoTheFit(List, BestSeeing, WorstSeeing);
BestImageConvolve();
return 1;
}



//! Carry out kernel fit.
/*!  Before entering here, 'this' should have the BestImage and
WorstImage fields assigned.  The various sizes of the stamps involved
in kernel fit are computed by OptParams::OptimizeSizes(). There is no
way yet to by-pass it and setup the fit conditions by yourself (it
would be simple to implement).  The returned value is 0 if the kernel
fit failed (and hence 'this' contains no viable kernel).  \param List
contains the stars elligible to fit the kernel.
*/


int KernelFit::DoTheFit(const BaseStarList &List, double &BestSeeing, double &WorstSeeing)
{
 optParams.OptimizeSizes(BestSeeing, WorstSeeing);
 optParams.OptimizeSpatialVariations(List.size());

 // sizes may have changed since last call so:
 if (BestImageStamps) delete BestImageStamps;
 BestImageStamps = new  
       StampList(*BestImage, List, optParams.HStampSize, optParams.MaxStamps);
 DeallocateConvolvedStamps();

 if (!BestImageStamps->size())
   {
     cerr << " No Stars to fit the kernel ... : giving up " << endl;
     KernAtCenterSum = 1;
     return 0;
   }

 cout << " starting kernel fit with 1/2 StampSize " << optParams.HStampSize 
      << " 1/2 KernelSize  " << optParams.HKernelSize 
      << " convolved Stamp size " << optParams.ConvolvedSize() << endl; 
 cout << " with " << BestImageStamps->size() << " stamps." << endl;
  clock_t tstart = clock();

  FitDifferentialBackground(3.0);

  KernelsFill();

  ComputeMAndB();
   cout << " finished computation of m and b" << endl;
  if (!Solve())
    {
    cout << " inversion failed  " << endl;
    delete BestImageStamps; BestImageStamps = NULL;
    DeallocateConvolvedStamps();
    return 0;
    }
  
  FilterStamps();
  clock_t tend = clock();
  cout << "CPU for the kernel fit " <<  float(tend- tstart)/float(CLOCKS_PER_SEC) << endl;
  DeallocateConvolvedStamps();
  if (BestImageStamps) { delete BestImageStamps; BestImageStamps = NULL;}
  return 1;
}


#ifdef USE_ROOT
/*
RUN_ROOTCINT
LINKDEF_CONTENT : #pragma link C++ class OptParams+;
LINKDEF_CONTENT : #pragma link C++ class XYPower+;
LINKDEF_CONTENT : #pragma link C++ class KernelFit+;

*/
#include "root_dict/kernelfitdict.cc"

#endif /* USE ROOT */
