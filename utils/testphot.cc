// -*- C++ -*-

#include <stdio.h> 
#include <getopt.h>

#include <iostream>
#include <cmath> // for PI

#define  CHECK_BOUNDS

#include "fitsimage.h"
#include "reducedimage.h"
#include "transformedimage.h"
#include "detection.h"
#include "apersestar.h"


extern "C" {volatile float rndm_();}

 #define rndm rndm_


static double sq(const double &x) { return x*x;}


static double GaussRndm()
{
  /*
  double r = 
  double phi = 2*M_PI*rndm();
  return r*cos(phi);
  */
  return sqrt(-2*log(rndm()))*cos(2*M_PI*rndm());
}


void usage()
{
  cerr << "usage: testphot [OPTIONS]" << endl;
  cerr << "where OPTIONS can be:" << endl;
  cerr << "   -N:   nimages   [default: 1]" << endl;
  exit(EXIT_FAILURE);
}



int main(int argc, char**argv)
{
  char c;
  int nImages = 1; // atoi(args[1]);
  int nx = 2000;
  int ny = 4000;
  
  while( (c=getopt(argc, argv, "N:h")) != -1 )
    switch(c) {
    case 'N':
      nImages = atoi(optarg);
      break;
    case 'h':
      usage();
      break;
    default:
      usage();
    }
  
  for (int k =0; k <nImages; ++k) {
    //    char name[16];
    stringstream sstrm;
    sstrm << "toc" << k;
    string name = sstrm.str();
    string namef = name + string("f");
    
    cout << " on attaque " << name << endl;
    
    remove(name.c_str());
    DbImage db(name.c_str());
    if (!db.IsValid()) db.Create("here");
    string elixirName = db.ElixirName();
    remove(elixirName.c_str());
    
    remove(namef.c_str());
    DbImage dbf(namef.c_str());
    if (!dbf.IsValid()) dbf.Create("here");
    string elixirNamef = dbf.ElixirName();
    remove(elixirNamef.c_str());
    
    
    FitsImage* I = new FitsImage(elixirName,nx,ny);
    FitsImage* If = new FitsImage(elixirNamef, nx, ny);
    double skysig = 10;
    (*I) += skysig*skysig;
    (*If) += skysig*skysig;
    
    for (int j=0; j<ny; j++)
      for (int i=0 ; i <nx; i++) {
	double v = skysig*GaussRndm();
    	(*I)(i,j) += v;
    	(*If)(i,j) += v;
      }
    
    double seeing=2;
    
    int step = 40;
    
    //  double rad[] = {2,3,4,5,6,7,8,9,10};
    // 3.23,5.67,7.12,9.76,11.56};
    //#define nrad (sizeof(rad)/sizeof(rad[0]))
    //      AperSEStarList toto;
    
    int stamp = 10;
    
    for (double x = step; x< nx -stamp; x+= step)
      {
	double fmax = 5*x/step;
	if (x+step > nx -stamp) // derniere colonne: saturation franche
	  {
	    fmax = max(65000., fmax*2);
	  }
	for (double y = step; y < ny - stamp ; y += step)
	  {
	    int imin = int(floor(x-stamp));
	    int imax = int(ceil(x+stamp));
	    int jmin = int(floor(y-stamp));
	    int jmax = int(ceil(y+stamp));
	    
	    for (int i=imin; i <=imax; ++i)
	      for (int j=jmin; j<=jmax; ++j)
		{
		  if (i<0 || i >=nx || j<0 || j >=jmax) continue;
		  double v = fmax*exp(-0.5*(sq(x-i) + sq(y-j))/sq(seeing));
		  //		  (*I)(i,j)+= v;
		  //		  (*If)(i,j)+= v;
		}
	  }
	
      } // fin boucle etoiles
    
    (*I)(500, 500) = 65000;
    (*If)(500, 500) = 65000;

    
    
    delete I;
    If->SetWriteAsFloat();
    delete If;
    
    string dbname(name);
    //    system(("elixir_to_toads "+dbname).c_str());
    //    system(("make_catalog -o -S "+dbname).c_str());
    //    system(("mkcat2 "+dbname).c_str());
    //    remove(db.FitsWeightName().c_str());
    //    remove(db.ElixirName().c_str());
    
    string dbnamef(namef);
    //    system(("elixir_to_toads "+dbnamef).c_str());
    //    system(("make_catalog -o -S "+dbnamef).c_str());
    //    system(("mkcat2 "+dbnamef).c_str());
    //    remove(dbf.FitsWeightName().c_str());
    //    remove(dbf.ElixirName().c_str());
    
    
  } // fin boucle images
  
  return EXIT_SUCCESS;
}
