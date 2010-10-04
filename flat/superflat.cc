#include <iostream>
//#include <fstream.h>
#include <cstdio>
#include <cmath>
#include <vector>

#include "fileutils.h"
#include "vutils.h"
#include "fitsimage.h"
#include "imageback.h"
#include "image.h"
#include "frame.h"
#include "superflat.h"
//#include "dbimage.h"
#include "fitsslice.h"
#include "fitsset.h"
#include "fitstoad.h"
#include "imageinterpolation.h"
#include "matvect.h"
#include "fringeutils.h"

#ifdef A_FAIRE

faire quelque chose de convenable pour le calcul de la saturation
actuellement, ca ne marche que pour le CFHT

render calling sequence of all routines the same kind 
for SkyFlat and Bias images (use pointers rather than references)
separate properly routines

Faire le menage : degager vieille routines, splitter superflat.cc (?).

#endif


#ifndef M_PI
#define  M_PI           3.14159265358979323846  /* pi */
#endif

  // uncomment this to use the new defringing method
#define NEW_DEFRINGING

static string calibdir(const string &InFileName)
{
char dirname1[256];
char dirname2[256];
DirName(InFileName.c_str(),dirname1); 
DirName(dirname1,dirname2);
string calib_dir = AddSlash(string(dirname2))+"calib/";
if (FileExists(calib_dir)) return calib_dir;
calib_dir =  AddSlash(string(dirname1))+"calib/";
if (FileExists(calib_dir)) return calib_dir;
return "";
}

//***************************************************************



//********************************************************
//! computes overscan and removes bias. Returns the average overscan value for all amps.
double BiasCorrect_and_trim(FitsImage &Current, const Image *Bias)
{
  Frame Illu = TotalIlluRegion(Current);
  // using Frame's for image sections: min=first pix, max = last pix
  // hence:
  int Nx_Image = int(Illu.Nx());
  int Ny_Image = int(Illu.Ny());
  
  if (Current.Nx()==Nx_Image && Current.Ny()==Ny_Image)
    {
      cout <<" Image "<< Current.FileName() 
	   <<" : bias, nlcorrections and trim -> already processed. Nothing done."
	   << endl;
      return 0.0;
    }

  double finalvalue = 0.0;
  int Namp;
  if (Current.HasKey("TOADNAMP")) Namp = Current.KeyVal("TOADNAMP");
  else Namp = 1;

  if (Current.HasKey("AMPNAME")) Namp = 1;

  int UseOnlyOverscan = ( !Bias );
 
  for (int i=1; i<= Namp; ++i)
    {
      int iamp = i;

      Frame overscanRegion =  OverscanRegion(Current, iamp);
      Pixel SigmaOverscan;
      Pixel MedianOverscan = Current.MedianInFrame(overscanRegion, 
						   SigmaOverscan);
      
      cout <<"Overscan region "<< iamp << " in image :   Median = " 
	   << MedianOverscan << "  Sigma = " 
	   << SigmaOverscan << endl;

      finalvalue += MedianOverscan;

      Frame illuRegion = IlluRegion(Current,iamp);
      int x0 = int(illuRegion.xMin); // belongs to image
      int y0 = int(illuRegion.yMin);
      int x1 = int(illuRegion.xMax);// beyond last
      int y1 = int(illuRegion.yMax);
     
      cout <<" x0 : "<< x0 << " y0 : " << y0 << " x1 : "<< x1 << " y1 : " << y1 << endl;

      
      if (UseOnlyOverscan == 0)
	{
	  Pixel SigmaBiasOverscan;
	  //correct image overscan with overscan in bias frame 
	  Pixel MedianBiasOverscan = Bias->MedianInFrame(overscanRegion, 
							 SigmaBiasOverscan);
	  cout <<"Overscan region "<< iamp << " in bias  :   Median = "
	       << MedianBiasOverscan << "  Sigma = " 
	       << SigmaBiasOverscan << endl;

	  const Image &bias = *(Bias); // handy 
	  Pixel correction = MedianOverscan - MedianBiasOverscan;
	  
	  for (int j=y0; j<y1; j++) for (int i=x0; i<x1; i++)
	    Current(i,j) = Current(i,j) - bias(i,j) - correction ;
  	}

      else
	  {
	    for (int j=y0; j<y1; j++) for (int i=x0; i<x1; i++)
	      Current(i,j) -= MedianOverscan;

	  }
    }
  Current.Trim(Illu);
  return finalvalue/Namp;
}


//***********************************************
//!  is a quick procedure to normalize an image
double ImageNormalisation(FitsImage &Current)
{
  int Namp;
  if (Current.HasKey("TOADNAMP")) Namp = Current.KeyVal("TOADNAMP");
  else Namp = 1;

  float sigma;
  float mean;

  if ((Namp==1) || !getenv("SEPARATE_AMPLI"))
    {
      Current.SkyLevel(&mean,&sigma);
      Current /= mean;
    }
  else
    {
      for (int i=1; i<= Namp; ++i)
	{
	  Frame illuRegion = AmpRegion(Current,i);
	  Current.SkyLevel(illuRegion,&mean,&sigma);
	  double norm = 1/mean;
	  Current.SubimageMultiply(illuRegion,norm);
	}
    }
  Current.SkyLevel(&mean,&sigma);
  cout << "Image normalised : mean = "<< mean <<"  sigma = "<< sigma << endl;
  return mean;
}

//************************************************
//! combines images to make a superflat
Image *MakeSuperFlat(const FitsSet &FitsFileSet, const Image *Bias, const Image *SkyFlat)
{
  vector<string> reducNames;
  
  //int nx_image = FitsFileSet.NxTot();
  //int ny_image = FitsFileSet.NyTot();
  int nx_image = FitsFileSet.Nx();
  int ny_image = FitsFileSet.Ny();
  int skyflat = (SkyFlat && SkyFlat->Nx() == nx_image  && SkyFlat->Ny() == ny_image);
  
  int nImages = FitsFileSet.size();
  if (nImages <2) 
    {
      cout <<" you requested to make a superflat from " << nImages << " images : we give up"  << endl;
      return NULL;
    }
  
  cout <<" WE USE A SIMPLE MEDIAN !!!!"<< endl;
  for(int k=0; k<nImages; k++)
    {
      string FitsName = FitsFileSet[k];
      string reducName = tempnam("/tmp",BaseName(FitsName).c_str());
      reducNames.push_back(reducName);
      if (!FileExists(reducName)) // stupid test 
	{
	  FitsImage *loop = new FitsImage(FitsName);
	  //	  FitsImage loop(FitsFileSet[k]);
	  FitsImage out(reducName,*loop,*loop);
	  delete loop;
	  BiasCorrect_and_trim(out, Bias);
	  if (skyflat) out /= (*SkyFlat);
	  else if(SkyFlat) 
	    {
	      cerr <<" The given skyflat doesn'y have the same nx or ny as the list images" << endl;
	    }
	  ImageNormalisation(out);
	}
    }
  
  Image *Flat = new Image(nx_image, ny_image);

  Image Count(nx_image, ny_image);

  int sliceSize = 100;
  
  int overlap = 0;
  FitsParallelSlices slices(sliceSize,overlap);
  
  for (int i=0; i<nImages; ++i) slices.AddFile(reducNames[i]);
  
  Pixel *pixelValues = new Pixel[nImages];
  
  cout << "Building a median master frame..." << endl;;
  do
    {
      //cout <<"slice : j = "<< slices.ImageJ(0) << endl; 
      for (int j=0; j<slices.SliceSize(); j++) for (int i=0; i<nx_image; i++)
 	{
	  for (int k=0; k<nImages; k++) pixelValues[k] = (*slices[k])(i,j);
	  int j_true = slices.ImageJ(j);
	  float mean;
	  int n = nImages;
	  mean = FArrayMedian(pixelValues, n);
	  (*Flat)(i,j_true) = mean;
	}
    }
  while (slices.LoadNextSlice());
  
  for (int i=0; i<nImages; ++i) remove(reducNames[i].c_str());
  delete [] pixelValues;
  
  if (skyflat) *Flat *= *SkyFlat;

  return Flat;
}

//****************************************************
//! produces a dead pixels mask image and smooths the flat with a box of 64
Image* DeadPixelImageAndFlatSmoothing(Image &Flat, double FlatMin, double FlatMax)
{
int Nx_Image = Flat.Nx();
int Ny_Image = Flat.Ny();
 cout << "Building a dead pixel image and smoothing the flat" << endl;
Image *DeadPic = new Image(Nx_Image,Ny_Image);

ImageBack back(Flat,64);

for (int j=0; j<Ny_Image; j++) for (int i=0; i<Nx_Image; i++)
  {
    double value = Flat(i,j);
    if (value<=FlatMin || value>=FlatMax)
      {
	(*DeadPic)(i,j)=1.0;

	if (value<=0.0)
	  {
	    Flat(i,j) = back.Value(i,j);

	    if (back.Value(i,j)<=0.0)
	      {
		Flat(i,j)=1.0;
		(*DeadPic)(i,j)=2.0;
	      }
	  }
      }
  }
return DeadPic;
}

//*******************************
//! returns the max value of an image and cut 
double ImageMaxCut(const Image &Current, const Image &Flat)
{
int nx=Current.Nx();
int ny=Current.Ny();

double maxValue=0.0;


/* if there is no flat specified, use only the central part of the image to calculate the max */
if (Flat.Nx()==0 || Flat.Ny()==0)
  {
    for (int j=ny/4; j<(3*ny/4); j++) for (int i=nx/4; i<(3*nx/4); i++)
      {
	double current=Current(i,j);
	if (maxValue<current) maxValue = current;
      }

    return maxValue;
  }

for (int j=0; j<ny; j++) for (int i=0; i<nx; i++)
  {
    double flatValue=Flat(i,j);
    double current=Current(i,j);
    if (flatValue<1.2 && flatValue>0.8 && maxValue<current)
      maxValue = current;
  }
for (int j=0; j<ny; j++) for (int i=0; i<nx; i++)
  {
    if (Current(i,j) > maxValue)
      Current(i,j) = maxValue;
  }
cout <<"MaxValue of Image is = "<< maxValue <<endl;
return maxValue;
}


static void CleanBadPix(FitsImage &Current, const Image &Flat, 
			double ReplacementVal)
{
  int namp = Current.KeyVal("TOADNAMP");
  for (int iamp = 1; iamp <= namp; ++iamp)
    {
      Frame frame = AmpRegion(Current,iamp);
      Pixel flatAverage, flatSig;
      Flat.SkyLevel(frame, &flatAverage, &flatSig);

      int starti = int(frame.xMin);
      int endi = int (frame.xMax);
      int startj = int(frame.yMin);
      int endj = int(frame.yMax);
      Pixel lowCut = 0.8*flatAverage;
      Pixel highCut = 1.2*flatAverage;
      for (int j=startj; j<endj; j++) for (int i=starti; i<endi; i++)
	{
	  Pixel flatValue=Flat(i,j);
	  if (flatValue<lowCut || flatValue>highCut) Current(i,j) = ReplacementVal;
	}
    } // end loop on amps
}


//*******************************
//! makes a raw average starting from an image list. The final image still has overscan (usefull for MasterBias making).
Image *MakeRawAverage(const FitsSet &FitsFileSet)
{

Image *Average = new Image(FitsFileSet.NxTot(), FitsFileSet.NyTot());
int nImages = FitsFileSet.size();

for(int i=0; i<nImages; i++)
  {    
    FitsImage loop(FitsFileSet[i]);
    *Average += loop;
  }

if (nImages != 0)
  {
    *Average /= nImages;
  }
 else cout << " MakeRawAverage : division by 0 !!!!!!! " << endl;

return Average;
}


Image *MakeRawAverageAndSigma(const FitsSet &FitsFileSet, Image *Sigma, const int normalize, const int normTime)
{
  //added the substraction of the mean overscan

  Image *Average = new Image(FitsFileSet.NxTot(), FitsFileSet.NyTot());
  int nImages = FitsFileSet.size();
  
  if (nImages <2) 
    {
      cout <<" you requested to make a master from " << nImages << " images : we give up"  << endl;
      return NULL;
    }
  
  int nx_image = FitsFileSet.NxTot();
  int sliceSize = 100;
  if (FitsFileSet.NyTot()<sliceSize) sliceSize = FitsFileSet.NyTot();

  int overlap = 0;
  FitsParallelSlices slices(sliceSize,overlap);
  double *norm = new double[nImages];
  Pixel *meanOS = new Pixel[nImages];
  Pixel *sigmaOS = new Pixel[nImages];

  for (int i=0; i<nImages; ++i)
    {
      slices.AddFile(FitsFileSet[i]);
      FitsImage Im(FitsFileSet[i]);
      int Namp=1;

      if (!Im.HasKey("WRITEDAT"))
	{
	  const Frame overscan= OverscanRegion(Im,Namp);   
	  meanOS[i] = Im.MedianInFrame(overscan, sigmaOS[i]);
	  cout <<" Will subtract mean overscan region for image "<< i <<" = "<< meanOS[i]<<endl;
	}
      else meanOS[i] = 0;

      if (normalize)
	{
	  if (normTime==1)
	    norm[i] = Im.KeyVal("TOADEXPO");
	  else 
	    {
	      float mean,sig;
	      Im.SkyLevel(&mean,&sig);
	      norm[i] = mean-meanOS[i];
	    }
	}
    }

  Pixel *pixelValues = new Pixel[nImages];
  
  cout << "Building a mean master frame from " << nImages <<" images" << endl;;
  do
    {
      for (int j=0; j<slices.SliceSize(); j++) for (int i=0; i<nx_image; i++) 
	{
	  for (int k=0; k<nImages; k++)
	    {
	      pixelValues[k] = (*slices[k])(i,j)-meanOS[k];
	      //	      if (normalize) pixelValues[k] *= norm[0]/norm[k];
	      if (normalize) pixelValues[k] /= norm[k];
	    }
	  int j_true = slices.ImageJ(j);
	  float mean;
	  float median;
	  float sigma;
	  Fmean_median_sigma(pixelValues, nImages, mean, median, sigma);
	  (*Average)(i,j_true) = mean;
	  (*Sigma)(i,j_true) = sigma;
	}
    }  
  while (slices.LoadNextSlice());
  
  delete [] pixelValues;
  delete [] norm;
  
  return Average;
}

//*******************************
//! makes a raw median starting from an image list. The final image still has overscan (usefull for MasterBias making).
Image *MakeRawMedian(const FitsSet &FitsFileSet, const int normTime)
{

  Image *Median = new Image(FitsFileSet.NxTot(), FitsFileSet.NyTot());
  int nImages = FitsFileSet.size();
  if (nImages <2) 
    {
      cout <<" you requested to make a master from " << nImages << " images : we give up"  << endl;
      return NULL;
    }
  
  int nx_image = FitsFileSet.NxTot();
  int sliceSize = 100;
  if (FitsFileSet.NyTot()<sliceSize) sliceSize = FitsFileSet.NyTot();

  int overlap = 0;
  FitsParallelSlices slices(sliceSize,overlap);
  double *timeExp = new double[nImages];
  for (int i=0; i<nImages; ++i)
    {
      slices.AddFile(FitsFileSet[i]);
      if (normTime)
	{
	  FitsHeader head(FitsFileSet[i]);
	  timeExp[i] = head.KeyVal("TOADEXPO");
	}
    }

  Pixel *pixelValues = new Pixel[nImages];
  
  cout << "Building a median master frame from " << nImages <<" images" << endl;;
  do
    {
      for (int j=0; j<slices.SliceSize(); j++) for (int i=0; i<nx_image; i++) 
	{
	  for (int k=0; k<nImages; k++)
	    {
	      pixelValues[k] = (*slices[k])(i,j);
	      if (normTime) pixelValues[k] /= timeExp[k];
	    }
	  int j_true = slices.ImageJ(j);
	  float mean;
	  mean = FArrayMedian(pixelValues, nImages);
	  (*Median)(i,j_true) = mean;
	}
    }  
  while (slices.LoadNextSlice());
  
  delete [] pixelValues;
  delete [] timeExp;
  
  return Median;
}

//*******************************
//! makes a raw masked median starting from an image list and an mask image list. The final image still has overscan (usefull for MasterBias making).
Image *MakeRawMaskedMedian(const FitsSet &FitsFileSet, const FitsSet &MaskFitsFileSet, const int normTime)
{
  int useMask=1;

  Image *Median = new Image(FitsFileSet.NxTot(), FitsFileSet.NyTot());
  int nImages = FitsFileSet.size();
  if (nImages <2)
    {
      cout <<" you requested to make a master from " << nImages << " images : we give up"  << endl;
      return NULL;
    }

  if (MaskFitsFileSet.size()!= FitsFileSet.size())
    {
      cout <<" not the same number of images and masks !!!"<<endl;
      return NULL;
    }

  int nx_image = FitsFileSet.NxTot();
  int sliceSize = 100;
  if (FitsFileSet.NyTot()<sliceSize) sliceSize = FitsFileSet.NyTot();

  int overlap = 0;
  FitsParallelSlices slices(sliceSize,overlap);
  FitsParallelSlices maskslices(sliceSize,overlap);
  double *timeExp = new double[nImages];
  for (int i=0; i<nImages; ++i)
    {
      slices.AddFile(FitsFileSet[i]);
      maskslices.AddFile(MaskFitsFileSet[i]);
      if (normTime)
	{
	  FitsHeader head(FitsFileSet[i]);
	  timeExp[i] = head.KeyVal("TOADEXPO");
	}
    }

  Pixel *pixelValues = new Pixel[nImages];
  
  cout << "Building a median master frame from " << nImages <<" images USING Mask Images "<< endl;;
  do
    {
      for (int j=0; j<slices.SliceSize(); j++) for (int i=0; i<nx_image; i++) 
	{
	  int npix = 0;
	  for (int k=0; k<nImages; k++)
	    {
	      if ((*maskslices[k])(i,j)>0)
		{
		  pixelValues[k] = (*slices[k])(i,j);
		  if (normTime) pixelValues[k] /= timeExp[k];
		  npix++;
		}
	    }
	  int j_true = slices.ImageJ(j);
	  float mean;
	  mean = FArrayMedian(pixelValues, npix);
	  (*Median)(i,j_true) = mean;
	}
    }  
  while (slices.LoadNextSlice());
  
  delete [] pixelValues;
  delete [] timeExp;
  
  return Median;
}

#include "histo1d.h"


#ifdef STORAGE
//*******************************************
/*! compute the most probable value of the pixel values (for the
 upper half of the distribution. It happens to be a bad 
"estimator" of the saturation 
*/

static double HistoEndMax(const Image &Current, double MaxValue)
{
  Histo1d histo(500,MaxValue/2., MaxValue);
  Pixel *p = Current.begin();
  Pixel *pend = Current.end(); 
  for ( ;p<pend; ++p) histo.Fill(*p);
  double satur;
  histo.MaxBin(satur);
  return satur;
}
#endif

//*****************************************
//! reads a gain file associated with an image
static int read_gains(const string InFileName, double *Gain, int Ccd /*daschan in fact */)
{
string gains = calibdir(InFileName) + "gains.dat";
FILE *file;
if (!(file = fopen(gains.c_str(),"r"))) return 0;

cout << "Found a file : " << gains << endl;
double gain = -1;
char line[256];
while (fgets(line,256,file))
  {
  int id;
  double gain_read;
  sscanf(line,"%i %lf",&id, &gain_read);
  if (id == Ccd) gain = gain_read;
  }
fclose(file);
if (gain != -1.)
  {
  *Gain = gain;
  return 1;
  }
else 
  {
    cerr << " could not find CCD # " << Ccd << " in file " << gains << endl;
    return 0;
  }
}
  
#include "alltelinst.h"

//*****************************************
//! sniffs if gain multiplied (by simply checking if gain is 1) and 
//! if not, multiplies it amp by amp
double GainMultiply(FitsImage &Current)
{
  int Namp;
  if (Current.HasKey("TOADNAMP")) Namp = Current.KeyVal("TOADNAMP");
  else Namp = 1;

  double gain = Current.KeyVal("TOADGAIN");
  if (gain==1)
    {
      cout << "Gain == 1 !!!! Image Already in photoelectrons ?"<< endl;
      return gain;
    }

  if ((Namp==1) || !getenv("SEPARATE_AMPLI"))
    {
      int chip = Current.KeyVal("TOADCHIP");
      if (!read_gains(Current.FileName(), &gain, chip))
	{
	  cout << "Converting image to photelectrons with " << gain << " e-/ADU" << endl;
	  Current *= gain;
	}
    }
  else
    {
      for (int i=1; i<= Namp; ++i)
	{
	  int chip = Current.KeyVal("TOADCHIP");
	  if (!read_gains(Current.FileName(), &gain, chip)) gain = AmpGain(Current,i);
	  if (gain != 1)
	    {
	      cout << "Converting Amp "<< i <<" to photelectrons with " << gain << " e-/ADU" << endl;
	      Frame illuRegion = AmpRegion(Current,i);
	      Current.SubimageMultiply(illuRegion,gain);
	    }
	}
    }
  Current.AddOrModKey("TOADGAIN",1.0," image (hopefully) in photoelectrons");
  Current.AddKey("OLDGAIN",gain," average gain used to convert pixels from adu to photoelectrons");
  double readnoise = Current.KeyVal("TOADRDON");
  Current.AddOrModKey("TOADRDON", readnoise*gain," original ro noise * gain");
  return gain;
}



// Modification du calcul de la saturation le 3 oct 2000 Julien
/* if anybody figures out how what this routine does, please write it here */
static double Saturation(const Image &image)
{
  float maxVal = image.MaxValue(); 
  double Saturation=maxVal, Sat;
  double scale = 50.0;
  int loop = 0;
  while (Saturation==maxVal && loop<2)
    {
      loop++;
      scale *= 10.0;
      Histo1d histo(int(maxVal/scale), 0, maxVal + 1 ); //HC
      for (int j = 0; j < image.Ny() ; ++j)
	for (int i = 0; i < image.Nx() ; ++i)
	  {
	    histo.Fill(image(i,j), 1 ); // HC
	  }
      
      double xMax;
      double maxContent = histo.MaxBin(xMax);
      double maxcoup = xMax;
      int count =1, bestMaxCount =0;
      for (int l=0; l<20; l++)
	{
	  histo.ZeroBins(maxVal*0.05*l, maxVal*0.05*(l+1));
	  Sat = xMax;
	  maxContent = histo.MaxBin(xMax);
	  if (Sat == xMax)
	    {
	      count ++;
	      if (count >= bestMaxCount && Sat != maxcoup)
		{ 
		  bestMaxCount = count;
		  Saturation =Sat;
		}
	    }
	  else 
	    {
	      count =1;
	    }
	}
    }
  return Saturation;
}



//*******************************************
//! does the flatfielding. 
/*! assumes that the flat file exists. If the bias does not exist, overscan 
is used. If the fringe exists, it is used for fringe subtraction.
File names are used rather than images, because we do not need all of them 
at the same time. Bug: if something goes wrong (different image sizes,
different chips or filters), the trimmed but unflatfielded image is written
as if everything was right. */
int FlatFieldImage(const string &InFileName, const string &FlatName, 
		   const string &BiasName, const string &FringeName, 
		   const string &OutFileName)
{

  // quick check on image and bias
  FitsImage *in = new FitsImage(InFileName);
  if (!in->IsValid())
    {
      delete in; return 0; // FitImage constructor produces a message.
    }
  
  FitsImage *bias = NULL;
  if (FileExists(BiasName)) 
    {
      bias = new FitsImage(BiasName);
      if (!bias->IsValid())
    {
      cerr << " FlatFieldImage : " << BiasName << " is not a valid FITS file \n";
      delete in; return 0;
    }
     if (!in->SameChip(*bias)) 
     { 
	 cerr << " FlatFieldImage : " << InFileName << " is not the same chip as " << BiasName << endl;
	 delete in; delete bias; return 0;
     }
     if (!in->SameSize(*bias))
	{
	  cout << "Bias frame "<< BiasName << " ain't got same size as " << InFileName << endl;
	  cout << "Will try with overscan " << endl;
	  delete bias;
	  bias = NULL;
	}
    }

  FitsImage outFits(OutFileName, *in, *in);
  // done with the input: 
  delete in; in = NULL; // to prevent accidental use.
  // switched to true if we reach the end of the routine:
  outFits.EnableWrite(false);

  double overscan = BiasCorrect_and_trim(outFits,bias);
  if (bias) delete bias; bias = NULL;// saves memory and prevent accidental use.
  
  if (getenv("WRITE_BIAS_REMOVED")) 
    // for debugging purposes : can save the image just after bias subtraction
    {
      string name = CutExtension(OutFileName)+"_br.fits";
      FitsImage(name, outFits, outFits);
    }

    // get the Image inside the FitsImage: ( not a copy !!) 
  Image &out = outFits;
  
  FitsImage  Flat(FlatName); 
  if (!Flat.IsValid())
    {
      cerr << " no flat " << FlatName << " found " << endl; 
      return 0;
    }
  if (!outFits.SameChipFilter(Flat))
    {
      cerr << " flat " << FlatName << " and image " << InFileName << endl 
	   << " do not refer to the same chip and/or filter " << endl;
      return 0;
    }

  // make sure the flat is trimed by triming it (why ?)
  // Flat.Trim(TotalIlluRegion(Flat));
  
  if (!outFits.SameSize(Flat))
    {
      cerr << " flat " << FlatName << " and image " << InFileName << endl 
	   << " do not have the same size " << endl;
      
      return 0;
    }
  
  // make sure the flat is normalized to 1 by renormalizing it
  // ImageNormalisation(Flat);
  // No : this has to be done somewhere else and once for all,
  // because there could be a good reason for flat averages different from 1,
  // such as different efficiencies of CCDs within a mosaic...
  
  // this is the flatfielding
  out /= Flat;
  
  if (getenv("WRITE_FLAT_FIELDED")) 
    // for debugging purposes : can save the image before fringe removal
    {
      string name = CutExtension(OutFileName)+"_ff.fits";
      FitsImage(name, outFits, out);
    }
  
  // remove a fringe pattern
  if (FileExists(FringeName))
    {
      float nsigma = 3;
      string band = StringToUpper(outFits.KeyVal("TOADBAND"));
      if(band == "Z")
	nsigma = 5;
      if(FringeUtils::RemoveFringes(outFits,FringeName,0,nsigma,0,1)!=0) {
	cout << "ERROR in superflat at FringeUtils::RemoveFringes" << endl;
      }
      outFits.AddOrModKey("FRINGED","SUBNEW","Fringe pattern subtracted (new method)");
    }
  
  // convert in photoelectrons
  /* No, because it is useless once one uses weights, and assumes that
   gains for different CCDs within a mosaic are correct
   (i.e. equalized by flatfielding). This is no the case for Megacam
   where the gains refer to the usual e/ADU, but flats compensate for
   both efficiency and actual electronic gain.  So a mosaic image
   nicely flatfielded for uniform response becomes non uniform if
   transformed to photelectrons, and for Megacam, these are NOT
   photoelectrons  */
  // double gain = GainMultiply(outFits);

  // computes max and saturation level
  Pixel sky,skysig;
  out.SkyLevel(&sky,&skysig);
  CleanBadPix(outFits,Flat,sky);
  double minValue = sky-7*skysig;
  //double satur = (65535-overscan)*gain;
  double satur = Saturation(outFits);
  
  cout << " enforce min and max at values : " << minValue << " and " << 1.05*satur << endl;
  out.EnforceMinMax(minValue,satur*1.05);  
  
  // fill in the header
  outFits.AddOrModKey("SATURLEV",satur," saturation value of the image after flat fielding");
  
  float sigmaflat=0.0,meanflat=0.0;
  Flat.MeanSigmaValue(&meanflat,&sigmaflat);
  outFits.AddKey("SIGFLAT",sigmaflat,"Sigma (in e-) of the flat used, gives an idea of saturation");
  
  float skylev, sigth=0.0, sigexp=0.0;
  out.MeanSigmaValue(&skylev,&sigexp);
  float readnoise = outFits.KeyVal("TOADRDON");
  if (skylev > 0.0) sigth = sqrt(skylev + readnoise*readnoise);
  outFits.AddKey("SKYLEV",skylev,"Sky Level in photoelectrons");
  outFits.AddKey("SKYSIGTH",sigth,"Sigma of Sky Level obtained from sqrt(sky+ronoise^2)");
  outFits.AddKey("SKYSIGEX",sigexp,"Sigma of Sky Level obtained from Image");
  
  outFits.AddCommentLine("flat fielding done using masterflat, Overscan removed");
  cout << "Final flatfielding for " << BaseName(outFits.FileName()) 
       << " : Sqrt(Sky+Readout^2) = " << sigth
       << "  RMS(image) = " << sigexp << endl;

  // flatfielding successful : write (i.e. keep) the result
  outFits.EnableWrite(true); 
  return 1;
}

/*This routine was cleaned up in 04/04, but not tested because we currently
  have no images that require this processing  */
void ImageAlreadyFlatFielded(const string &InFileName, const Image &Flat, 
                     const string &OutFileName)
{
FitsImage in(InFileName);
if (!in.IsValid()) return;

FitsImage outImage(OutFileName, in, in);


// converts in photoelectrons, if it's already in photoelectrons, nothing done
// No : this prooves to be a bad idea
// double gain = GainMultiply(outImage);

Image &out = outImage;
 
double maxValue = ImageMaxCut(out,Flat);
double satur =  Saturation(out);

outImage.AddKey("MAXPIX",maxValue,"max value of the image after flat fielding");
outImage.AddKey("SATURLEV",satur,"max value of the image after flat fielding");

double sigmaflat = 0.02; /* it's an average value, as there is no flat specified */
outImage.AddKey("SIGMAFLAT",sigmaflat,"Default value put by hand: no flat");

Pixel skylev=0.,sigexp=0.0;
out.MeanSigmaValue(&skylev,&sigexp);
// correct even if image is GainMultiplied:
 double readnoise = outImage.KeyVal("TOADRDON"); 
double gain = outImage.KeyVal("TOADGAIN");

 double sigth=0.0;
 if (skylev > 0.0) sigth = sqrt(skylev/gain + readnoise*readnoise);
outImage.AddKey("SKYLEV",skylev,"Sky Level in image units");
outImage.AddKey("SKYSIGTH",sigth,"Sigma of Sky Level obtained from sqrt(sky+ronoise^2)");
outImage.AddKey("SKYSIGEX",sigexp,"Sigma of Sky Level obtained from Image");

outImage.AddCommentLine("Already flatfielded without TOADS, Overscan removed");
}

#ifdef STORAGE

void  FilterFringes(Image &FringePattern, const int GridSize, const int NInc)
{
  // TO OPTIMIZE GREATLY
  int nx = FringePattern.Nx();
  int ny = FringePattern.Ny();
  int halfsize = (GridSize-1)/2; 
  double *colarr = new double[GridSize]; // the column array 
  double mean,disp,medcol,dispmax; //stats
  double dMaxInc = M_PI/NInc; // maximum inclination
  double bestangle = 0; 
  double x,y;
  int row,col,nrow,ncol; // scripts for the grid
  int mins = -halfsize/2-1;
  int maxs = halfsize/2+1;
  
  cout << "Sweeping out a grid of " << GridSize << " X " << halfsize << endl;
  //loop on pixels
  for (int j=0; j<ny; j++) for (int i=0; i<nx; i++)
    {
      cout << "\r" << i << "  " << j << "  before = " << FringePattern(i,j) << flush ; 
      dispmax = 0;
      //loop on angles
      for (double angle=0;angle<M_PI;angle+=dMaxInc)
	{
	  mean = 0;
	  disp = 0;
	  ncol = 0;
	  // loop on grid columns
	  for (col=mins;col<=maxs;col++,ncol++)
	    {
	      nrow = 0;
	      //loop on rows
	      for (row=-halfsize;row<=halfsize;row++,nrow++)
		{
		  x = i + row*cos(angle) - col*sin(angle);
		  y = j + row*sin(angle) + col*cos(angle);
		  if ((x < 1) || (x > nx) ||  (y < 1) || (y > ny) ) 
		    {
		      colarr[nrow] = FringePattern(i,j);
		    }
		  else 
		    {
		      colarr[nrow] = FringePattern(int(x+0.5),int(y+0.5));
		    }
		}
	      medcol = DArrayMedian(colarr,GridSize);
	      disp += medcol*medcol ;
	      mean += medcol; 		
	    }
	  mean /= ncol;
	  disp = disp/ncol-mean*mean;
	  disp = (disp > 0) ? disp : 0;
	  if (disp>dispmax) {dispmax=disp;bestangle = angle;}
	}      	      
      nrow = 0;
      // filter pixel i,j by rotating with found angle 
      for (row=-halfsize;row<=halfsize;row++,nrow++)
	{
	  x = i + row*cos(bestangle);
	  y = j + row*sin(bestangle);
	  if ((x < 1) || (x > nx) ||  (y < 1) || (y > ny) ) 
	    {
	      colarr[nrow] = FringePattern(i,j);
	    }	
	  colarr[nrow] = Interpolate(FringePattern,x,y);
	}
      
      FringePattern(i,j) = DArrayMedian(colarr,GridSize);
      cout << "  after  = " << FringePattern(i,j) << flush;
    }
  cout << "done!!" << endl;
  delete [] colarr;
	  
}

void  RawHighFreqFilter(Image &Img, const int GridSize)
{
  int nx = Img.Nx();
  int ny = Img.Ny();
  Image temp = Img;

  for (int j=0; j<ny; j++) for (int i=0; i<nx; i++)
    {
      double sum = 0;
      int npix = 0;
      for (int k=0; k<GridSize; k++) for (int l=0; l<GridSize; l++)
	{
	  int u = i+k-(GridSize/2);
	  int v = j+l-(GridSize/2);
	  if (u<0 || u>nx || v<0 || v>ny) continue;
	  sum += Img(u,v);
	  npix++;
	}
      temp(i,j) = sum/npix;
    }
  Img = temp;
}


Image* MakeFringePattern(const Image &FringedImage, const Image &BlankImage,const double &Nsig, const bool &BackFlag, const bool &FilterFlag)
{
  int ny = FringedImage.Ny(), nx = FringedImage.Nx();
  Image *Fringe = new Image(nx,ny);

   *Fringe = FringedImage / BlankImage;
   
   float skylev,skysigma;
   Fringe->SkyLevel(&skylev,&skysigma);
   //   Fringe->MeanSigmaValue(&skylev,&skysigma);
   Fringe->EnforceMinMax(skylev-Nsig*skysigma,skylev+Nsig*skysigma);

   *Fringe -= 1.0;
   if (BackFlag) 
       {
         int meshStep = 256;
         cout << " Removing background of mesh step " << meshStep <<endl;
         ImageBack back(*Fringe,meshStep);
         for (int j=0; j<ny; j++) for (int i=0; i<nx; i++)
            {
                (*Fringe)(i,j) -= back.Value(i,j);
            }
       }
   if (FilterFlag)
       {
          int GridSize = 5;
          cout << "Applying low pass band filter of grid size "<< GridSize << endl;
          RawHighFreqFilter(*Fringe,GridSize); 
       }
   return Fringe;
}


static double NormalizeFringes( const Image &TheImage, const Image& Fringes)
{
if (!TheImage.SameSize(Fringes))
  {
    cerr << "NormalizeFringes :  images do not have the same size " << endl;
    return 0;
  } 
float sky, sig_sky;
TheImage.MeanSigmaValue(sky,sig_sky);

}

static double NormFactor(const Image &Im1,const Image &Im2, const double NSig)
  // Returns the empirical covariance of 2 images.
{
  if (Im1.Nx() != Im2.Nx() || Im1.Ny() != Im2.Ny() ) 
    {
      cerr << "NormFactor : Image 1 and 2 are not the same size. Returning 0" << endl;
      return 0;
    }
  float mean1,sig1,mean2,sig2;
  Im1.MeanSigmaValue(&mean1, &sig1);
  Im2.MeanSigmaValue(&mean2, &sig2);
  
  double sum1 = 0;
  double sum2 = 0;
  double sum12 = 0;
  double sum22 = 0;
  int npix = 0;
  Pixel *p1, *p2;
  Pixel *end = Im1.end();
  double cut1  = NSig*sig1;
  double cut2  = NSig*sig2;

  int jump = 10;
  for (p1=Im1.begin(), p2=Im2.begin(); p1<end; p1 += jump, p2 += jump)
    {
    if (fabs(*p1-mean1) > cut1) continue;  
    if (fabs(*p2-mean2) > cut2) continue;  
      sum1 += *p1;
      sum2 += *p2;
      sum12 += (*p1)*(*p2);
      sum22 += (*p2)*(*p2);
      npix++;
    }
  sum1 /= npix;
  sum2 /= npix;
  sum12 /= npix;
  sum22 /= npix;
  return ((sum12-sum1*sum2)/(sum22 - sum2*sum2)); 
}


#include "kernelfit.h" /* for XYPower */



#include "vutils.h" /* for MatSolve */

static double NormFactor(const Image &Im,const Image &F, const double NSig,
const int BackDegree)
  /* returns the fringes coefficient, computed as the value that mimizes 
     the sky variance, after subtraction of a poynomial of degree Deg in X,Y
     from the image (which is however NOT changed).
  */
{
  
  if (Im.Nx() != F.Nx() || Im.Ny() != F.Ny() ) 
    {
      cerr << "NormFactor : Image 1 and 2 are not the same size. Returning 0" << endl;
      return 0;
    }
  XYPower backModel(BackDegree);
  int nterms = backModel.Nterms();
  int msize = nterms+1;
  Mat A(msize,msize);
  Vect B(msize);
  Vect monom(nterms);

  Pixel imMean,imSig,fMean,fSig;
  Im.SkyLevel(&imMean, &imSig);
  F.SkyLevel(&fMean, &fSig);
  
  double cut1  = NSig*imSig;
  double cut2  = NSig*fSig;

  int jump = 10;
  for (int j=0 ; j< Im.Ny() ; j+= jump)
  for (int i=0 ; i< Im.Nx() ; i+= jump)
    {
      Pixel p1 = Im(i,j);
      if (fabs(p1-imMean) > cut1) continue;  
      Pixel p2 = F(i,j);
      if (fabs(p2-fMean) > cut2) continue;
      for (int q1=0; q1<nterms; ++q1)
	{
	  monom(q1) = backModel.Value(double(i), double(j),q1);
	

  for (int q2 = q1; q2<nterms; ++q2)
	    A(q1,q2) += monom(q1)*monom(q2);
	  B(q1) += monom(q1)*p1;
	  A(q1,nterms) += monom(q1)*p2;
	}
      B(nterms) += p1*p2;
      A(nterms,nterms)  += p2*p2;
    }
  /* symetrize */
  for (int q1=0; q1<nterms; ++q1) 
    for (int q2 = q1+1; q2<nterms; ++q2) A(q2,q1) = A(q1,q2); 
  if (!MatSolve(&A(0,0),msize,&B(0)))
    {
      cerr << " could not compute fringes normalization: no fringe subtraction" << endl;
      return 0;
    }
  return B(nterms);
}

double SurfaceFit(const Image &Im, const int BackDegree)
{
  XYPower backModel(BackDegree);
  int nterms = backModel.Nterms();
  int msize = nterms;
  //  int msize = nterms+1;
  Mat A(msize,msize);
  Vect B(msize);
  Vect monom(nterms);

  Pixel imMean,imSig;
  Im.SkyLevel(&imMean, &imSig);
  
  for (int j=0 ; j< Im.Ny() ; j++)
    for (int i=0 ; i< Im.Nx() ; i++)
      {
	Pixel p1 = Im(i,j);
	for (int q1=0; q1<nterms; ++q1)
	  {
	    monom(q1) = backModel.Value(double(i), double(j),q1);
	    
	    for (int q2 = q1; q2<nterms; ++q2)
	      A(q1,q2) += monom(q1)*monom(q2);
	    B(q1) += monom(q1)*p1;
	    //	    A(q1,nterms) += monom(q1)*p2;
	  }
	//	B(nterms) += p1*p2;
	//	A(nterms,nterms)  += p2*p2;
      }
  /* symetrize */
  for (int q1=0; q1<nterms; ++q1) 
    for (int q2 = q1+1; q2<nterms; ++q2) A(q2,q1) = A(q1,q2); 
  if (!MatSolve(&A(0,0),msize,&B(0)))
    {
      cerr << " could not compute image surface fit !!!!" << endl;
      return 0;
    }
  return B(nterms);
}


void RemoveFringes(Image &FringedImage, const Image &FringeMap,const bool &NormFlag)
{
  if (FringeMap.Nx() != 0 && FringeMap.Ny() !=0 )
  {
    cout << " Removing fringes" << endl;
    if (NormFlag)
      {
	double ncut = 3.0;
	double norm = NormFactor(FringedImage,FringeMap,ncut,3);
	cout << " Scale fringe factor = " << norm << endl;
        FringedImage -= norm*FringeMap; 
      }
    FringedImage -= FringeMap; 
  }
}

#endif
