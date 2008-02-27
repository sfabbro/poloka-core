#include <iostream>

#include "fitsimage.h"
#include "hbook.h"

#define NWPAWC 100000
float pawc_[NWPAWC];

#define MAXVAR 200
#define MAX_LENGTH 30
#define TOPDIR "//TOPDIR"

int open_hbook_file(const string &name)
{
int istat=0;
char toto2[256];
sprintf(toto2,name.c_str());
 int lrec = 1024; 
HROPEN(1,"TOPDIR",toto2,"N",lrec,istat);
if (istat == 0) return 1;
return 0;
}


static void usage(const char * pg)
{
  printf(" syntax : %s <image1> <image2> [<hbookfile>] [-n <number of pixels>\n",pg);
  exit(1);
}

int main(int argc, char **argv)
{
  vector<string> fileNames;
  int nPixTot = 100000;
  for (int i=1; i<argc;++i)
    {
      char * arg = argv[i];
      if (arg[0] != '-')
	{
	  fileNames.push_back(argv[i]);
	}
      else
	switch (arg[1])
	  {
	  case 'n' : ++i; nPixTot = atoi(argv[i]); break;
	  default : usage(argv[0]);
	  }
    }// end loop on args
  


  if (fileNames.size() !=2 && fileNames.size() != 3) usage(argv[0]);
  string nomhb = (fileNames.size() == 2) ? "truc.hbk" : fileNames[2];
  string nomim1 = fileNames[0];
  string nomim2 = fileNames[1];
  FitsImage img1(nomim1);
  FitsImage img2(nomim2);

  if ( ( img1.Nx() != img2.Nx() ) || ( img1.Ny() != img2.Ny() ) )
    {
      cerr << " Images de taille differente ! " << endl ;
      return 0 ;
    }

  HLIMIT(NWPAWC);
  int Id = 1 ;
  int dim = 4;
  char **tags = (char**) calloc(MAXVAR,sizeof(char*));
  tags[0] = "x" ;
  tags[1] = "y" ;
  tags[2] = "f1" ;
  tags[3] = "f2" ;
  //tags[4] = "fond" ;
  //tags[5] = "floc" ;
  //tags[6] = "flux" ;
  //tags[7] = "area" ;
 
  if (!open_hbook_file(nomhb))
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


 // A CHANGER AVEC DIM !!!!
   int npix =  img1.Nx()*img1.Ny();
   int nx = img1.Nx();
   int pas = (int) (npix/(1.*nPixTot));
   if( nPixTot > npix )
     {
       nPixTot = npix ; 
       pas = 1 ;
     }
  cout << "pas " << pas << endl ;
  float x[8];
  register Pixel *p;
  p = img1.begin() ;
  Pixel *start = p;
  int count =0;
  for ( int i=0;i<nPixTot;p += pas, i++)
    {
      int jpix=(p-start)/nx;
      int ipix=(p-start)%nx;
      x[0]=ipix;
      x[1]=jpix;
      x[2]=img1(ipix,jpix);
      x[3]=img2(ipix,jpix);
      HFN(Id,x);
      count++;
    }


 int icycle = 0;
 char toto5[50]=" ";
 HROUT(0,icycle,toto5);
 HREND(toptop);

 cout << "written " << count << " pixels " << " to " << nomhb << endl;


 return 0;
}

  
