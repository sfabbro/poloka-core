#include <iostream>
#include <string>
#include "fitsimage.h"
#include "hbook.h"
#include "fileutils.h"

#define NWPAWC 100000
float pawc_[NWPAWC];

#define MAXVAR 200
#define MAX_LENGTH 30
#define TOPDIR "//TOPDIR"

int open_hbook_file(string name)
{
  int istat=0;
  char toto1[50]="TOPDIR";
  char toto2[256];
  sprintf(toto2,name.c_str());
  char toto3[50]="N";
  int lrec= 1024;
  HROPEN(1,toto1,toto2,toto3,lrec,istat);
  if (istat == 0) return 1;
  return 0;
}

static void usage(const std::string &Prog)
{
  cerr << " usage : " << std::endl 
	   << Prog << " [-n npixtot] <image> [<hbbokfile>]"
       << std::endl;
  exit(1);
}


int
main(int argc, char **argv)
{
  std::string nomhb;
  std::string  imageName;
  int nPixTot = 50000;
  for (int i=1; i< argc; ++i)
    {
      char * arg = argv[i];
      if (arg[0] != '-')
	{
	  if (imageName == "") { imageName = argv[i]; continue;}
	  else {nomhb = argv[i]; continue;}
	}
      else
	switch (arg[1])
	  {
	  case 'n' : ++i; nPixTot = atoi(argv[i]); break;
	  default : usage(argv[0]);
	  }
    }
  if (imageName == "") usage(argv[0]);
  if (nomhb == "") nomhb = SubstituteExtension(imageName, ".hbk");

  FitsImage img(imageName);
  HLIMIT(NWPAWC);
  int Id = 1 ;
  int dim = 3;
  char **tags = (char**)calloc(MAXVAR,sizeof(char*));
  tags[0] = "x" ;
  tags[1] = "y" ;
  tags[2] = "f" ;
 
  if (!open_hbook_file(nomhb.c_str()))
    {
      printf( " could not open %s\n",nomhb.c_str()); return 0;
    }


  char ttags[MAXVAR][MAX_LENGTH];
  for (int i=0; i<dim; i++) 
    {strncpy(ttags[i],tags[i],MAX_LENGTH-1);}
  char title[128];
  strcpy(title,"toto");  
  char toptop[50]="TOPDIR";
  HBOOKN(Id,title,dim,toptop,60000,ttags);

  float x[3];
  /*
    int pas = 10 ;
    for(int i = 0 ; i < img.Nx() ; i+=pas)   
    for(int j = 0 ; j < img.Ny() ; j+=pas)
    {
    x[0]=i;
    x[1]=j;
    x[2]=img(i,j);
    HFN(Id,x);
    }*/

  int npix =  img.Nx()*img.Ny();
  int nx = img.Nx();
  int pas = (int) (npix/(1.*nPixTot));
  if( nPixTot > npix )
    {
      nPixTot = npix ; 
      pas = 1 ;
    }
  cout << "pas " << pas << endl ;

  register Pixel *p;
  p = img.begin() ;
  Pixel *start = p;
  for ( int i=0;i<nPixTot;p += pas, i++)
    {
      int jpix=(p-start)/nx;
      int ipix=(p-start)%nx;
      x[0]=ipix; 
      x[1]=jpix; 
      x[2]= *p;
      HFN(Id,x);
    }

  int icycle = 0;
  char toto5[50]=" ";
  HROUT(0,icycle,toto5);
  HREND(toptop);


  return 0;
}

  
