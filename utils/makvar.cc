#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <fstream>
#include <string>

#include "sub.h"

#include "fileutils.h"

#include "fitsimage.h"

#include <getopt.h>

#define MAXVAR 400

char* split_line(char *Line, float *X, const int Dim, int *nread)
{
  const char *p1;
  char* p2;
  int i;
  *nread = 0;
  if (strlen(Line) <= 1) return Line;
  memset(X,0,Dim*sizeof(X[0]));
  p1 = Line;
  for (i=0; i< Dim; ++i)
    {
    float value = strtod(p1,&p2);
    if (p2 == p1) break;
    X[i] = value;
    p1 = p2;
    }
  return p2;
};



int factoriel(int N)
{
  if(N<=1)
    return 1;
  else
    return(N*factoriel(N-1));
}


double hermite(int n, double x) {

  double h;

  switch (n) {
  case 0 :
    h = 1.;
      break;
  case 1 :
    h = 2*x;
      break;
  case 2 :
    h = 4.*x*x-2.;
      break;
  case 3 :
    h = 8*x*x*x-12.*x;
      break;
  case 4 :
    h = 16*pow(x,4)-48.*x*x+12.;
      break;
  case 5 :
    h = 32*pow(x,5)-160.*x*x*x+120.*x;
      break;
  case 6 :
    h = 64.*pow(x,6)-480.*pow(x,4)+720.*x*x-120.;
      break;
  case 7 :
    h = -1680*x + 3360*pow(x,3) - 1344*pow(x,5) + 
      128*pow(x,7);
      break;
  case 8 :
    h = 1680 - 13440*pow(x,2) + 13440*pow(x,4) - 
      3584*pow(x,6) + 256*pow(x,8);
      break;
  case 9 :
    h = 30240*x - 80640*pow(x,3) + 48384*pow(x,5) - 
      9216*pow(x,7) + 512*pow(x,9);
      break;
  case 10 :
    h = -30240 + 302400*pow(x,2) - 403200*pow(x,4) + 
      161280*pow(x,6) - 23040*pow(x,8) + 1024*pow(x,10); 
      break;
  case 11 :
    h = -665280*x + 2217600*pow(x,3) - 1774080*pow(x,5) + 
      506880*pow(x,7) - 56320*pow(x,9) + 2048*pow(x,11);
      break;
  case 12 :
    h = 665280 - 7983360*pow(x,2) + 13305600*pow(x,4) - 
           7096320*pow(x,6) + 1520640*pow(x,8) - 
           135168*pow(x,10) + 4096*pow(x,12);
      break;
  case 13 :
    h = 17297280*x - 69189120*pow(x,3) + 69189120*pow(x,5) - 
           26357760*pow(x,7) + 4392960*pow(x,9) - 
           319488*pow(x,11) + 8192*pow(x,13);
      break;
  case 14 :
    h = -17297280 + 242161920*pow(x,2) - 484323840*pow(x,4) + 
           322882560*pow(x,6) - 92252160*pow(x,8) + 
           12300288*pow(x,10) - 745472*pow(x,12) + 
           16384*pow(x,14);
      break;
#if 0
  case 15 :
    h = -518918400*x + 2421619200*pow(x,3) - 
   2905943040*pow(x,5) + 1383782400*pow(x,7) - 
   307507200*pow(x,9) + 33546240*pow(x,11) - 
      1720320*pow(x,13) + 32768*pow(x,15);
      break;
  case 16 :
    h = 518918400 - 8302694400*pow(x,2) + 
   19372953600*pow(x,4) - 15498362880*pow(x,6) + 
   5535129600*pow(x,8) - 984023040*pow(x,10) + 
   89456640*pow(x,12) - 3932160*pow(x,14) + 
      65536*pow(x,16);
      break;
  case 17 :
    h = 17643225600*x - 94097203200*pow(x,3) + 
   131736084480*pow(x,5) - 75277762560*pow(x,7) + 
   20910489600*pow(x,9) - 3041525760*pow(x,11) + 
   233963520*pow(x,13) - 8912896*pow(x,15) + 
      131072*pow(x,17);
      break;
  case 18 :
    h = -17643225600 + 317578060800*pow(x,2) - 
   846874828800*pow(x,4) + 790416506880*pow(x,6) - 
   338749931520*pow(x,8) + 75277762560*pow(x,10) - 
   9124577280*pow(x,12) + 601620480*pow(x,14) - 
      20054016*pow(x,16) + 262144*pow(x,18);
      break;
  case 19 :
    h = -670442572800*x + 4022655436800*pow(x,3) - 
   6436248698880*pow(x,5) + 4290832465920*pow(x,7) - 
   1430277488640*pow(x,9) + 260050452480*pow(x,11) - 
   26671841280*pow(x,13) + 1524105216*pow(x,15) - 
      44826624*pow(x,17) + 524288*pow(x,19);
      break;
  case 20 :
    h = 670442572800 - 13408851456000*pow(x,2) + 
   40226554368000*pow(x,4) - 42908324659200*pow(x,6) + 
   21454162329600*pow(x,8) - 5721109954560*pow(x,10) + 
   866834841600*pow(x,12) - 76205260800*pow(x,14) + 
   3810263040*pow(x,16) - 99614720*pow(x,18) + 
      1048576*pow(x,20); 
      break;
#endif


      default :
      h = 0;

  }

  return h;

}

double phin(int n1, int n2, double x, double y) {

  double pi = 3.1415926;
  double phi = hermite(n1,x)*hermite(n2,y)*exp(-0.5*(x*x+y*y))/
               sqrt(pow(2.,(n1+n2))*sqrt(pi)*factoriel(n1)*factoriel(n2));

  return phi;
};



void ComputeShapecoef(double xpix, double ypix, int xc, int yc , Image &Outwindow, const Image &InImage, double coef[20][20],int nmax,double seseeing)
{

  //cout << " xc Yc " << xc << ' ' << yc << endl;

  int width = 20;
  int minx = xc - width/2;
  int miny = yc - width/2;

  double back =0;  

  for (int n1=0;n1<nmax;n1++) {

    for (int n2=0;n2<nmax;n2++) {

      coef[n1][n2]=0; 
  

      for (int j=0 ; j< width; ++j)
	{
	  int jin = miny+j;
	  if (jin<0 || jin >= InImage.Ny()) continue; 
	  for (int i=0 ; i< width; ++i)
	    {
	      int iin = minx+i;
	      if (iin<0 || iin >= InImage.Nx()) continue;

	      double f = InImage(iin,jin)-back;
	      double beta = seseeing*1.1;

              //double x = (iin-xc)/beta;
              //double y = (jin-yc)/beta;

	    double x = (iin-xpix)/beta;
	    double y = (jin-ypix)/beta;

              //cout << f << "---" << phin(n1,n2,x,y) << endl;
	      coef[n1][n2] += f*phin(n1,n2,x,y)/beta;

	      //Outwindow(i,j) = InImage(iin,jin)-back;
	      //Outwindow(i,j) = phin(0,0,x,y);
	    }
	}
      //cout << "coefs (" <<n1 << "," << n2 << ") = " <<  coef[n1][n2]/coef[0][0] << endl;

    }
  }


}




void ComputeBadalign(int xc, int yc , int width, const Image &InImage, double back, double &v1, double &v2, double &v3)
{

  //cout << " xc Yc " << xc << ' ' << yc << endl;

  int minx = xc - width/2;
  int miny = yc - width/2;

  
  double r2in=0;
  double r5in=0;
  double r8in=0;
  double r11in=0;
  double rwin=0;

  for (int j=0 ; j< width; ++j)
    {
      int jin = miny+j;
      if (jin<0 || jin >= InImage.Ny()) continue; 
      for (int i=0 ; i< width; ++i)
	{
	  int iin = minx+i;
	  if (iin<0 || iin >= InImage.Nx()) continue;
	  rwin += InImage(iin,jin)-back;
	  if (sqrt(pow((double)(xc-iin),2)+pow((double)(yc-jin),2)) < 11) {
	    r11in += InImage(iin,jin)-back;
	    if (sqrt(pow((double)(xc-iin),2)+pow((double)(yc-jin),2)) < 8) {
	      r8in += InImage(iin,jin)-back;
	      if (sqrt(pow((double)(xc-iin),2)+pow((double)(yc-jin),2)) < 5) {
		r5in += InImage(iin,jin)-back;
		if (sqrt(pow((double)(xc-iin),2)+pow((double)(yc-jin),2)) < 2) {
		  r2in += InImage(iin,jin)-back;
		}
	      }
	    }
	  }
	  //Outwindow(i,j) = InImage(iin,jin)-back;
	}
    }

  //cout << r2in << " " << (r5in-r2in)/r2in << " " << (r8in-r5in)/r2in << " " 
  //     << (r11in-r8in)/r2in << " " << r11in/rwin << endl;

  v1= (r5in-r2in)/r2in;
  v2= (r8in-r5in)/r2in;
  v3= (r11in-r8in)/r2in;

}


void ComputeVars(string inName, int xcand,int ycand, 
		 double &v1,double &v2,double &v3) {

 
  int outSize = 40;
 
  FitsImage inImage(inName);
  //FitsHeader h(inImage);

  //double back = h.KeyVal("BACKSUM");
  //  h.ModKey("NAXIS1",outSize);
  // h.ModKey("NAXIS2",outSize);

  //Image vignetteImage(outSize, outSize); 


  ComputeBadalign(xcand,ycand,outSize,inImage,0,v1,v2,v3);


}




int
main(int argc, char **argv)
{

  string Pathname = argv[1];
  char*  Sseseeing = argv[2];
double seseeing = atof(Sseseeing);


  string DetFilename  = Pathname + "/det.list";
  string ImgFilename  = Pathname + "/calibrated.fz";
  string Det2Filename  = Pathname + "/det2.list";

  string outName = "stamp.fits";


 float x[MAXVAR];
 char line[8192];
 int count = 0;
 int miss = 0; int more = 0;

 int Dim = 32;
 if (Dim>MAXVAR)
   {
     Dim = MAXVAR;
   }

 int nmax = 4;
 double coef[20][20];
 
 FitsImage inImage(ImgFilename);

 int outSize = 20;
 Image vignetteImage(outSize, outSize); 

 //FitsHeader h(inImage);
  //double back = h.KeyVal("BACKSUM");

 //h.ModKey("NAXIS1",outSize);
 //h.ModKey("NAXIS2",outSize);


 FILE *in;
 FILE *out;

 char *vname;
 char **tags;
 int dim;
 dim=33;
 in = fopen(DetFilename.c_str(),"r");
 out = fopen(Det2Filename.c_str(),"w");

 if (!in) 
  {
  printf(" cannot open %s\n",DetFilename.c_str());
  }
 cout << DetFilename.c_str() << endl;
 if (!out) 
  {
  printf(" cannot open %s\n",Det2Filename.c_str());
  }
 cout << Det2Filename.c_str() << endl;
 

while (fgets(line,8192,in))
  {
  char *left_over;
  int nread;

  if (strlen(line) <= 1) {
    fputs(line,out);
    continue;
  }

  if (line[0] == '#') {
    if (strncmp(line,"# end",5) == 0) {
      fputs("# Disc var 1\n",out);
      fputs("# Disc var 2\n",out);
      fputs("# Disc var 3\n",out);
    }
    fputs(line,out);
    continue;
  }


  left_over = split_line(line,x,Dim,&nread);
  if (nread < Dim) miss++;
  if (atof(left_over)) more++;
  
  double xpix=x[0];
  double ypix=x[1];

  //  cout << " db1 " << endl;
  int xc = (int)xpix;
  int yc = (int)ypix;

  double v1;
  double v2;
  double v3;
  //ComputeVars(ImgFilename,xc,yc,v1,v2,v3);

  ComputeShapecoef(xpix,ypix,xc,yc,vignetteImage,inImage,coef,nmax,seseeing);

    // faire ici loop nmax * nmax pour output coef[i][j]
 
  char buffer[20];
  line[strlen(line)-1]=0;

  for (int i=0;i<nmax;i++) {
    for (int j=0;j<nmax;j++) {
                                                                               
      memset (buffer,0,20);
      gcvt (coef[i][j],10,buffer);
      buffer[strlen(buffer)] = ' ';
      buffer[strlen(buffer)] = 0;
      strcat(line,buffer);

    }
  }

  memset (buffer,0,20);
  buffer[strlen(buffer)] = '\n';
  buffer[strlen(buffer)] = 0;
  strcat(line,buffer);
                                                                                
  fputs(line,out);

  count++;
  //cout << count << endl;
  }

 fclose(out);

 //FitsImage outImage(outName,h,vignetteImage);

}











