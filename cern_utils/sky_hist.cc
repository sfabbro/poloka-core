
#include <iostream>
#include <stream.h>
#include <fstream.h>
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <limits.h>
#include <string.h>
#include <string>

#include <time.h>

#include "hbook.h"

#include "fitsimage.h"

#define NWPAWC 100000
float pawc_[NWPAWC];

//#define NOTDEAD

static FILE* kumac;

static char* BareFileName(char *FileName)
{
  char *p = FileName+strlen(FileName) - 1;
  while (p>=FileName && *p != '/') p--;
  return p+1;
}



int image_histos(char *FileName, char *Mask)
{
  if (!FileName) return 0;



  FitsImage img(FileName);
  static int id =11;
  static int idmask =101;

  int nx = img.Nx();
  int ny = img.Ny();
  int ntot = nx*ny;
  //double gain=img.KeyVal("GAIN");

  FitsImage *pmask = NULL;
  if ( Mask != NULL && (strlen(Mask) > 1 ) )        
    {
      pmask = new FitsImage(Mask) ;
      if ( pmask->Nx() != nx ||  pmask->Ny() != ny )
	{
	  cerr << "taille mask != taille image " << endl ;
	  return 0;
	}
    }





  char histName[80]; /* with HBOOK, have to pass a real array, not a pointer as far as I remember */
  char histMName[80]; /* with HBOOK, have to pass a real array, not a pointer as far as I remember */

  //  strcpy(histName,BareFileName(FileName));
  strcpy(histName,FileName);
  string temp1  = histName;
  string temp = "mask" + temp1 ;
  strcpy(histMName, temp.c_str());
  //Image *image = &img;
  //*image *= gain;

  double average = 0;
  double variance = 0;
  double value;
  double sigma = 0; 
  float aa, ss ;
  img.SkyLevel(&aa,&ss);
  average = aa ;
  sigma = ss ;
  cout << "image "<< histName  << " ave " 
       << average << " sqrt(ave) " << sqrt(average) 
       << " sigma " << sigma << " count " << ntot << endl;

  //  if (average > 2.0)
  char *nloo = getenv("NLOOP");
  int Nloop = 3 ;
  if ( nloo != NULL)
    {
      sscanf(nloo,"%d", &Nloop);
      cout << "N loop : " << Nloop << endl ;
    }
    
  int Nsig = 3.5 ;

  if (true)
    {
      for(int loop=0; loop< Nloop; loop++)
	{
	  int count=0;
	  double mean = average;
	  average = 0.0;
	  variance = 0.0;
	  for (int j=0; j<ny; j++) for (int i=0; i<nx; i++) 
	    {
	      value = img(i,j);
	      if (fabs(value-mean)<Nsig*sigma)
		{
		  average += value;
		  variance += value*value;
		  count++;
		}
	    }
	  average /= count;
	  variance = (variance/(count-1)) - (average*average);
	  sigma = sqrt(variance);
	  cout << "image "<< histName  << " ave " << average 
	       << " sqrt(ave) " << sqrt(average) << " sigma " 
	       << sigma << " count " << count << endl;
	}

      HBOOK1(id,histName,50, average-5.*sigma, average+5.*sigma, 0);
      if (pmask != NULL )
	HBOOK1(idmask,histMName,50, average-10.*sigma, average+10.*sigma, 0);
    }
  else
    {
      HBOOK1(id,histName,50, -0.5, 1.5, 0);
      if (pmask != NULL )
	HBOOK1(idmask,histMName,50, -0.5, 1.5, 0);
    }

  for (int i=0; i<nx; i++) for (int j=0; j<ny; j++)
    {
      HFILL(id,img(i,j),0.,1.);
      //  if (img(i,j)<=0.1)
      //  cout << "img("<<i<<","<<j<<") = "<<img(i,j)<<endl;
      if (pmask != NULL )
	if ((*pmask)(i,j)< 1.e-10)
	  HFILL(idmask,img(i,j),0.,1.);

    }


  fprintf(kumac,"h/plot %d\n",id);
  id++;
  if (pmask != NULL )
    {
      fprintf(kumac,"h/plot %d\n",idmask);
      idmask++;
    }
  return 1;
}


int image_histos(char *FileName)
{
  char *Mask ;
  return image_histos(FileName, Mask) ;

}






int main(int nargs,char **args)
{
 char *mask = getenv("MASK");
 



  string nom;

  int istat;
  HLIMIT(NWPAWC);
#define HBK_FILE_NAME "histos.hbk"

  char toto1[50]="TOPDIR";
  char toto2[50]="histos.hbk";
  char toto3[50]="N";

  HROPEN(1,toto1,toto2,toto3,1024,istat);
  if (istat !=0)
    {
      cout << "cannot open histos.hbk" << endl;
      exit(2);
    }

  kumac = fopen("plot_all.kumac","w");
  if (!kumac) 
    {
      cout << " cannot open plot_all.kumac ..." << endl;
      exit (2);
    }
  fprintf(kumac,"h/file 1 %s\n",HBK_FILE_NAME);
  int i;
  for (i=1 ; i<nargs; i++)
    {
      cout << " processing image " << i <<" : "<< args[i] << " --- histo : " << i+10 <<endl;
      image_histos(args[i], mask);
    }

  int icycle;
  char toto4[50]=" ";
  HROUT(0,icycle,toto4);
  HREND(toto1);
  

  return 1;

}
