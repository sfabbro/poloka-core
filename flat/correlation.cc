#include <iostream>
#include <string>
//#include <fstream>

#include "fitsimage.h"
//#include "frame.h"
#include "fitstoad.h"
//#include "alltelinst.h"
#include "vutils.h"

static int greatest_common_divider(int a, int b)
{
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

static void usage(const char *prog)
{
  cout << prog << " #npix <fitsfile(s)> " << endl; 
  exit(EXIT_FAILURE);
}  

int main(int argc,char **args)
{
  if (argc <3) {
    usage(args[0]);
  }
  int correlationlength = atoi(args[1]);
 
  for (int i=2; i< argc; ++i) {
    string imagename = args[i];
    FitsImage image(imagename);
    
    if (!image.IsValid()) {
      cerr << " skylev : invalid file : "  << args[i] << endl;
      continue;
    }
    
    Pixel mean_tot,sigma_tot;
    image.SkyLevel(&mean_tot,&sigma_tot);
    float cutoff = sigma_tot*4;
    
    int nx = image.Nx();
    int ny = image.Ny();
    int npix =nx*ny;
    int nvalues=100000;
    int step = npix/nvalues;
    
    if (npix < nvalues) {
      step = 1;
      nvalues = npix;
    }
    while (greatest_common_divider(nx,step) != 1)
      step--;
    nvalues = (npix/step);
    double *pixel = new double[nvalues]; // P(i)
    double *pixelprod = new double[nvalues*2]; // P(i)*P(i+correlationlength)
    
    int count=0;
    int count_prod=0;
    
    int i=0;
    int j=0;
    float val1,val2;
    
    while(true)
      {
	val1=image(i,j);
	if(fabs(val1-mean_tot)<cutoff) {
	  pixel[count++]=val1;
	  if(j+correlationlength<ny) {
	    pixelprod[count_prod++] = val1*image(i,j+correlationlength);
	  }
	  if(i+correlationlength<nx) {
	      pixelprod[count_prod++] = val1*image(i+correlationlength,j);
	  }
	}
	int toto = i+step;
	j += toto/nx;
	i = toto%nx;
	if (j>=ny-correlationlength) break;
      }
    
    double rms;
    double mean = clipmean(pixel, count, rms);
    double rms_cov;
    double mean_prod = clipmean(pixelprod, count_prod, rms_cov);
    double covariance = mean_prod-mean*mean;
    //cout << mean << " " << mean_prod << " " << covariance << " " << rms << endl;
    cout << imagename << " " << covariance/(rms*rms) << endl;
    /*
    ofstream stream("toto.list");
    stream << "# prod :" << endl;
    stream << "# i :" << endl;
    stream << "# end" << endl;
    for(int i=0;i<count_prod;i++) {
      stream << pixelprod[i] << " " << i << endl;
    }
    stream.close();
    */
    delete [] pixel;
    delete [] pixelprod;
    
  }
  return EXIT_SUCCESS;

}
