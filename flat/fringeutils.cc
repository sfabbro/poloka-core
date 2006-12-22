/* 
 * $Source: /cvs/snovae/toads/poloka/flat/fringeutils.cc,v $
 * $Revision: 1.2 $
 * $Author: guy $
 * $Date: 2006/12/22 13:35:40 $
 * $Name:  $
 */

#include <iostream>

// TOADS
#include "fitsimage.h"
#include "fitsimagearray.h"
#include "vutils.h"
#include "fringeutils.h"
#include "imageback.h"
#include "matvect.h"

//uncomment this to get a lot of output
#define VERBOSE 

// get cvs version of the code
#define CVSVERSION "$Revision: 1.2 $"


int FringeUtils::GreatestCommonDivider(int a, int b) {
  if (a<b) swap(a,b);
  do {
    int q = a/b;
    int r = a-b*q;
    if (r==0) return b;
    a = b;
    b = r;
  } while ( b!= 1);
  return 1;
}

bool  FringeUtils::IsSameSize(const Image &Img1, const Image &Img2) {
  return (Img1.Nx() == Img2.Nx()) && ( Img1.Ny() == Img2.Ny());
}



double FringeUtils::ScalarProduct(Image &Img1, Image &Img2, float nsigma, const Image* deadimage) {
  if(!IsSameSize(Img1,Img2)) {
    cout << "ERROR fringefinder::ScalarProduct Images do not have the same size" << endl;
    return 0;
  }
  if(deadimage)
    if(!IsSameSize(Img1,*deadimage)) {
      cout << "ERROR fringefinder::ScalarProduct Image and deadimage do not have the same size" << endl;
      return 0;
    }
  
  Pixel mean1,sigma1,max1;
  Pixel mean2,sigma2,max2; 
  Pixel val1,val2;

  
  Img1.SkyLevel(&mean1,&sigma1);
  Img2.SkyLevel(&mean2,&sigma2);
  
  max1=nsigma*sigma1;
  max2=nsigma*sigma2;
  
  
  Pixel *pix1 = Img1.begin();
  Pixel *pix1end = Img1.end();
  Pixel *pix2 = Img2.begin();
  Pixel *deadpix = 0;
  if(deadimage)
    deadpix = deadimage->begin();
  
  double sum = 0;
  unsigned int npix = 0;
  if(!deadimage) {
    if(nsigma <= 0 ) {
      npix = Img1.Nx()*Img1.Ny();
      for (; pix1 < pix1end; ++pix1, ++pix2) {
	sum += (*pix1-mean1)*(*pix2-mean2);
      }
    } else { 
      for (; pix1 < pix1end; ++pix1, ++pix2) {
	val1 = *pix1-mean1;
	val2 = *pix2-mean2;
	if( (fabs(val1)<max1) && (fabs(val2)<max2) ) {
	  sum += val1*val2;
	  npix ++;
	}
      }
    }
  }else{
    if(nsigma <= 0 ) {
      for (; pix1 < pix1end; ++pix1, ++pix2, ++deadpix) {
	if( (*deadpix==0)) {
	  sum += (*pix1-mean1)*(*pix2-mean2);
	  npix ++;
	}
      }
    }else{
      for (; pix1 < pix1end; ++pix1, ++pix2, ++deadpix) {
	val1 = *pix1-mean1;
	val2 = *pix2-mean2;
	if( (*deadpix==0) && (fabs(val1)<max1) && (fabs(val2)<max2) ) {
	  sum += val1*val2;
	  npix ++;
	}
      } 
    }
  }
  
#ifdef VERBOSE 
  if(npix)
    cout << sum/npix << " ((" << npix << ") fraction of pixels removed = " << 100.*(1-((float)npix)/(Img1.Nx()*Img1.Ny())) << " %)" << endl; 
  else
    cout << "no valid pixel" << endl;
#endif
  
  if(npix)
    return (sum/npix);
  else
    return 0;
} 

void  FringeUtils::ClearImage(Image &Img, const double NSigma) {
  
  bool keepbelow = true;
  double nsigma = NSigma;
  
  if(nsigma<0) {
    keepbelow = false;
    nsigma = - nsigma;
    cout << "!!! Warning , nsigma <0, so I will keep signal ABOVE fabs(nsigma) !!!" << endl;
  }
  float mean,sigma;
  Img.SkyLevel(&mean,&sigma);
  
  if(sigma<=0)
    return;
  
  int nx = Img.Nx();
  int ny = Img.Ny();
  for (int j=0; j<ny; j++)
    for (int i=0; i<nx; i++)
      if(fabs((Img(i,j)-mean)/sigma)>nsigma) {
	if(keepbelow)
	  Img(i,j)=0;
      }else{
	if(!keepbelow)
	  Img(i,j)=0;
      }
}


void  FringeUtils::SmoothFilter(Image &Img) {
  int nx = Img.Nx();
  int ny = Img.Ny();
  Image temp = Img;
  double sum;
  int npix;
  int x,y;
  for (int j=0; j<ny; j++) {
    for (int i=0; i<nx; i++) {
      sum = 0;
      npix = 0;
      for (int k=-1; k<2; k++) {
	x = i+k;
	if(x<0 || x>=nx)
	  continue;
	for (int l=-1; l<2; l++) {
	  y = j+l;
	  if(y<0 || y>=ny)
	    continue;
	  npix++; // neighbouring pixel is in the image
	  
	  if(x==0 && y==0) {
	    sum += 0.5*Img(x,y);
	    continue;
	  }
	  if(x==0 || y==0) { 
	    sum += 0.1*Img(x,y);
	    continue;
	  }
	  sum += 0.025*Img(x,y);
	}
      }
      if(npix)
	temp(i,j) = sum/npix;
    }
  }
  Img = temp;
}

Mat  FringeUtils::ScalarProductMatrix(vector<string> &filelist, const Image* deadimage) {
  
  int nImages = filelist.size();
  if (nImages <2) {
    cout <<" you requested to make a ScalarProductMatrix with less than 2 images !!" << endl;
    return Mat(0,0);
  }
  
  Mat scalar_product_matrix(nImages,nImages);
  
  double scalarproduct;
  //string fitsnamei,fitsnamej;
  //FitsImage *imagei,*imagej;

#ifdef VERBOSE  
  int ntot = nImages*(nImages+1)/2;
#endif 
  float nsigma = 3;
  
  int n=0;
  // need to optimize memory allocation ??
  for(int i = 0 ;i<nImages;i++) {
    FitsImage imagei(filelist[i].c_str(),RO);
    for(int j = i ;j<nImages;j++) {
      FitsImage imagej(filelist[j].c_str(),RO);
      n++;
#ifdef VERBOSE
      cout << "(" << i << "," << j << ") ("<< n << "/" << ntot << ") => " ;
#endif 
      scalarproduct = ScalarProduct(imagei,imagej,nsigma,deadimage);
      scalar_product_matrix(i,j)=scalarproduct;
      scalar_product_matrix(j,i)=scalarproduct;
    }
  }
  return scalar_product_matrix;
}

string FringeUtils::GetVersion() {
  string version = CVSVERSION;
  return version;
}

void FringeUtils::Normalize(Image &img,const float nsigma_cut) {
  float mean,sigma;
  img.SkyLevel(&mean,&sigma);

#ifdef VERBOSE
  cout << "before: mean = " << mean << ", sigma = " << sigma << endl;
#endif 
  
  for(int i=0;i<img.Nx();i++)
    for(int j=0;j<img.Ny();j++)
      img(i,j)=img(i,j)-mean;
  
  double scalarp = ScalarProduct(img,img,nsigma_cut);
  if(scalarp<=0) {
    cout << "ERROR in FringeUtils::Normalize : ScalarProduct = " << scalarp << endl;
    return;
  }
  for(int i=0;i<img.Nx();i++)
    for(int j=0;j<img.Ny();j++)
      img(i,j)=img(i,j)/sqrt(scalarp);

#ifdef VERBOSE
  img.SkyLevel(&mean,&sigma);
  cout << "after: mean = " << mean << ", sigma = " << sigma << endl;
#endif
}

void FringeUtils::ClippedAverageRms(double *pixelValues, const int count,
				double &average, double &rms) {
  
  double median = DArrayMedian(pixelValues, count);
  
  average = median;
  rms = 100000.0;
  
  double meanval,sigval;
  int nloop = 10;
  double quantization = 1.;
  for(int loop=0; loop<nloop; loop++)
    {
      meanval = 0.0;
      sigval = 0.0;
      int nval = 0;
      for (int i=0; i < count; i++)
	{
	  double val = pixelValues[i] - average;
	  if (i!=0) 
	    // compute difference to previous value (they are sorted) to evaluate data quantization.
	    {
	      double diff = pixelValues[i] - pixelValues[i-1];
	      if (diff != 0) quantization = min(diff,quantization);
	    }
	  if (fabs(val)<4*rms) 
	    { // Fix for images with very small sigma (eventually 0)
	      meanval += val;
	      sigval += val*val;
	      nval++;
	    }
	}
      /* meanval is the difference to average, see just above */
      meanval /= nval;
      average += meanval;
      sigval = (sigval/(nval-1)) - (meanval*meanval);
      if (sigval > 0)
	{
	  sigval = sqrt(sigval);
	  rms = sigval;
	}
      else sigval = sqrt(sigval/(nval-1));
      if (meanval == 0) break; // why?
      if (rms < quantization) break; // no need to loop
    }
}

int FringeUtils::RemoveFringes(FitsImage &image, const string &fringefilename, int nvec , float nsig, bool substractbg, bool verbose) {
  
  if(!image.IsValid()) {
    cout << "ERROR in FringeUtils::RemoveFringes : Input image is not valid" << endl; 
    return -1;
  }
  if(image.FileMode()!=RW) {
    cout << "ERROR in FringeUtils::RemoveFringes : Input image must be writable" << endl; 
    return -1;
  }
  if(verbose)
    cout << " Opening fringes file " << fringefilename << endl;
  
  FitsImageArray fringe_map(fringefilename,RO);
  if(!fringe_map.IsValid()) {
    cout << "ERROR in FringeUtils::RemoveFringes : " << fringefilename << " is not valid" << endl;
    return -1;
  }
  
  
  Pixel sky,sigma;
  if(verbose) {
    image.SkyLevel(&sky, &sigma);
    cout << " Before removing fringes : Sky = " << sky << "  Sigma = " << sigma << endl;
  }
  ImageBack *back = 0;
  if(substractbg) { 
    if(verbose)
      cout << " Creating background image ..." << endl;
    back = new ImageBack(image,500); 
    if(verbose)
      cout << " Removing background ..." << endl;
    image -= *(back->BackgroundImage());
    if(verbose) {
      image.SkyLevel(&sky, &sigma); 
      cout << " Without background : Sky = " << sky << "  Sigma = " << sigma << endl;
    }
  }
  double norme=0;
  
  
  
  int nvectot = 1;
  nvectot = fringe_map.NHDU();
  //if(fringe_map.HasKey("NEXTEND"))
  //nvectot = (int)fringe_map.KeyVal("NEXTEND");
  
  if(verbose)
    cout << " Number of fringes vectors in this file = " << nvectot << endl;
  
  if(nvec==0)
    nvec=nvectot;
  if(nvec>nvectot)
    nvec = nvectot;
  
  char comment[256];
  image.AddCommentLine("----------------------------------------");
  sprintf(comment,"Processed with defringe2 version %s",CVSVERSION);
  image.AddCommentLine(comment);
  
  sprintf(comment,"Finding fringes with %d vectors",nvec);
  image.AddCommentLine(comment);

  sprintf(comment,"Fringe file : %s",fringefilename.c_str());
  image.AddCommentLine(comment);
  
  if(verbose)
    cout << " Finding fringes with " << nvec << " vectors" << endl;
  
  FitsImageArray::Status status;
  for(int ivec=0; ivec<nvec;ivec++) {
    status = fringe_map.At(ivec+1);
    if(status!=FitsImageArray::OK) {
      cout << " Problem reading vector " << ivec << endl;
      continue;
      //break;
    }
    
    if(verbose)
      cout << ivec << ": " << fringe_map.KeyVal("EXTVER") << " " << fringe_map.KeyVal("EXTNAME") << " ..." << endl;
    
    norme = FringeUtils::ScalarProduct(image,fringe_map,nsig)/FringeUtils::ScalarProduct(fringe_map,fringe_map,nsig);
    if(verbose)
      cout << " coef = " << norme;

    sprintf(comment,"Fringe coefficient for vector %d = %f",ivec,norme);
    image.AddCommentLine(comment);
     
    image -= norme*fringe_map;
    image.SkyLevel(&sky, &sigma); 
    cout << " (=> Sky = " << sky << "  Sigma = " << sigma << ")" << endl;
  }
  
  image.AddCommentLine("----------------------------------------");
  
  if(substractbg && back) {
    cout << " Adding background ..." << endl;
    image += *(back->BackgroundImage());
    delete back;
  }

  image.AddOrModKey("SKYLEV",sky,"Sky Level in photoelectrons");
  image.AddOrModKey("SKYSIGEX",sigma,"Sigma of Sky Level obtained from Image");
  image.AddOrModKey("FRINGED","SUBNEW","Fringe pattern subtracted (second method)");
  
  return 0;
}

 
bool FringeUtils::IsDefringed(const FitsImage &image) {
  
  if(image.HasKey("IMCMB_FT"))
    return true;
  if( !image.HasKey("FRINGED"))
    return false;
  
  string defringe_status = image.KeyVal("FRINGED");
  return defringe_status=="SUB" || defringe_status=="SUBNEW";
}
  

bool FringeUtils::IsDefringedWithNewMethod(const FitsImage &image) {
  
  if( !image.HasKey("FRINGED"))
    return false;
  
  string defringe_status = image.KeyVal("FRINGED");
  return defringe_status=="SUBNEW";
}
  

bool FringeUtils::IsANewFringePattern(const string &filename) {
  FitsImage image(filename,RO);
  if( ! image.HasKey("EXTNAME") )
    return false;
  
  string extname = image.KeyVal("EXTNAME");
  return  ((int)extname.find("fringe")) >= 0 ;
}
