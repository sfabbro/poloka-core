#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <iomanip.h>
#include <fstream.h>
#include <getopt.h>
#include <string.h>
#include <vector.h>

#include "dodetection.h"
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


int
main(int argc, char **argv)
{
 string nomhb;
 if (argc == 3)
   {
     nomhb = "truc.hb";
   }
 if (argc == 4) nomhb = argv[3];
 if (argc < 3 || argc > 4)
   {
     printf(" syntax : im2tup <image1> <image2> <hbbokfile> \n");
     exit(1);
   }
 string nomim1 = argv[1];
 string nomim2 = argv[2];
 FitsImage img1(nomim1);
 FitsImage img2(nomim2);

 if ( ( img1.Nx() != img2.Nx() ) || ( img1.Ny() != img2.Ny() ) )
   {
     cerr << " Images de taille differente ! " << endl ;
     return 0 ;
   }

 HLIMIT(NWPAWC);
 int Id = 1 ;
 int dim = 8;
 char **tags = (char**) calloc(MAXVAR,sizeof(char*));
 tags[0] = "x" ;
 tags[1] = "y" ;
 tags[2] = "f1" ;
 tags[3] = "f2" ;
 tags[4] = "fond" ;
 tags[5] = "floc" ;
 tags[6] = "flux" ;
 tags[7] = "area" ;
 
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

 double rad_flux = 2.5* 1.7 ;

 // A CHANGER AVEC DIM !!!!
 float x[8];
 for(int i = 0 ; i < img1.Nx() ; i+=10)   
   for(int j = 0 ; j < img1.Ny() ; j+=10)
     {
       x[0]=i;
       x[1]=j;
       double rad1 = 6.*rad_flux/2.5 ; 
       double rad2 = 8.*rad_flux/2.5 ; 
       //double ffd = FondLocal( img1 , i, j, rad1, rad2, 0.15, 0.15);
       int area;
       double fluxloc = Flux_Aperture(img1, i, j,rad_flux ,0,area );
       double flux0 = Flux_Aperture(img1, i, j,rad_flux ,0.,area );
       x[2]=img1(i,j);
       x[3]=img2(i,j);
       x[4]= 0;
       x[5]= fluxloc;
       //x[6]= flux0;
       x[6]= 0;
       x[7]= area;
       HFN(Id,x);
     }


 int icycle = 0;
 char toto5[50]=" ";
 HROUT(0,icycle,toto5);
 HREND(toptop);


 return 0;
}

  
