#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <fstream.h>
#include <getopt.h>
#include <string.h>

#include "fileutils.h"
#include "fitsimage.h"
#include "image.h"
#include "dbimage.h"

int
Taille_Scale(int Taille, int delta);
Image *
Rescale(Image & image, int delta);
void 
Add(Image const & image1, Image & imagetot, int xor_tot, int yor_tot) ;
Image *
Retourne(Image *pim1);
Image *
Retourne5(Image *pim1);

 void usage();
void usage()
{
  cerr << "usage: " << endl ; 
  cerr << "-n : binning (default = 4)" << endl ;
  cerr << "-t: tag d'un des ccd sans le _numero" << endl ;
  cerr << "-T: tag d'un des ccd (pas forcement le 1 ), version < 31/08/99" << endl ;
  cerr << "OU BIEN : " << endl ;
  cerr << "-1: nom image ccd1 " << endl ;
  cerr << "-2: nom image ccd2 " << endl ;
  cerr << "-3: nom image ccd3 " << endl ;
  cerr << "-4: nom image ccd4 " << endl ; 

  cerr << "OU BIEN : " << endl ;
  cerr << "-i: nom image 4xfits, sans le .fits, version > 31/08/99, recherche alors les imahes nom_1.fit, nom_2.fit ...." << endl ;

  cerr << endl
       << "on les place ainsi  DASCHAN (TOADSCHIP): " << endl 
       << endl 
       << " o__________ " << endl
       << " |         |" << endl
       << " |   3     |" << endl
       << " -----------" << endl
       << " o---------- -----o" << endl
       << " |         | |    |" << endl
       << " |   4     | |    |" << endl
       << " ----------- |5(2)|" << endl
       << " o---------- |    |" << endl
       << " |         | |    |" << endl
       << " |   1     | |    |" << endl
       << " ----------- ------" << endl<< endl;
  exit(-1);
}
int
main(int argc, char **argv)
{
  char c;
  string nom1, nom2, nom3, nom4, nomtag , nomim;
  int ntot = 0 ;
  bool ouitag = false , ouioldtag = false, ouiim = false;
  int delta = 4 ;

  while ((c = getopt(argc, argv, "h1:2:3:4:t:T:i:n:")) != -1) 
    {
      switch (c)
	{
	case 'h' :
	  usage();
	  break;
	  case 'n' :
	  delta = atoi(optarg) ;
	  break;
	  case '1' :
	  nom1 = optarg ; ntot++;
	  break;
	  case 'i' :
	   nomim = optarg ;
	   ouiim = true ;
	  break;
	  case '2' :
	  nom2 = optarg ;ntot++;
	  break;
	  case '3' :
	  nom3 = optarg ;ntot++;
	  break;
	  case '4' :
	  nom4 = optarg ;ntot++;
	  break;
	  case 't' :
	  nomtag = optarg ;
          ouitag = true ;
	  break;
	  case 'T' :
	  nomtag = optarg ;
          ouitag = true ;
          ouioldtag = true ;
	  break;
	default:
	  usage();
	}
    }
  if (optind >= argc+1) usage();

  /* on les place ainsi  DASCHAN (TOADSCHIP):

        o__________
        |         |
        |   3     |
        -----------
        o---------- -----o
        |         | |    |
        |   4     | |    |
        ----------- |5(2)|
        o---------- |    |
        |         | |    |
        |   1     | |    |
        ----------- ------

*/



  if ( ouitag && ouioldtag )
    {
      int num = 0 ;
      sscanf(nomtag.c_str(), "r%d", &num) ;
      cout << "num " << num << endl ;
      
      DbImage dbim(nomtag.c_str());
      FitsHeader header(dbim.FitsImageName(Calibrated));
      int ccd = header.KeyVal("DASCHAN");
      if ( ccd == 5 )
        num -= 1 ;
      if ( ccd == 3 )
        num -= 2 ;
      if ( ccd == 4 )
        num -= 3 ;

      char c[100] ;
      sprintf(c, "r%d", num) ;

      DbImage dbim1(c);
      nom1 =  dbim1.FitsImageName(Calibrated) ;

      num += 1 ;
      sprintf(c, "r%d", num) ;
      DbImage dbim2(c);
      nom2 = dbim2.FitsImageName(Calibrated) ; 
 
      num += 1 ;
      sprintf(c, "r%d", num) ;
      DbImage dbim3(c);
      nom3 = dbim3.FitsImageName(Calibrated) ; 

      num += 1 ;
      sprintf(c, "r%d", num) ;
      DbImage dbim4(c);
      nom4 = dbim4.FitsImageName(Calibrated) ; 

      cout << nom1 << endl
	   << nom2 << endl
	   << nom3 << endl
	   << nom4 << endl ;
    }

  if ( ouitag && !(ouioldtag) )
    {
      char c[256] ;
      sprintf(c, "_%d", 1) ;
      string s = nomtag + c ;
      DbImage dbim1(s);
      nom1 =  dbim1.FitsImageName(Calibrated) ;

      
      sprintf(c, "_%d", 2) ;
      s = nomtag + c ;
      DbImage dbim2(s);
      nom2 = dbim2.FitsImageName(Calibrated) ; 
 
      sprintf(c, "_%d", 3) ;
      s = nomtag + c ;
      DbImage dbim3(s);
      nom3 = dbim3.FitsImageName(Calibrated) ; 

      sprintf(c, "_%d", 4) ;
      s = nomtag + c ;
      DbImage dbim4(s);
      nom4 = dbim4.FitsImageName(Calibrated) ; 

      cout << nom1 << endl
	   << nom2 << endl
	   << nom3 << endl
	   << nom4 << endl ;
    }

 
  if ( ouiim )
    {
      char c[256] ;
      sprintf(c, "_%d", 1) ;
      nom1 = nomim + c + ".fit" ;
 
      sprintf(c, "_%d",2 ) ;
      nom2 = nomim + c + ".fit" ; 
 
      sprintf(c, "_%d",3 ) ;
      nom3 = nomim + c + ".fit" ; 
 
      sprintf(c, "_%d",4 ) ;
      nom4 = nomim + c + ".fit" ; 

      cout << nom1 << endl
	   << nom2 << endl
	   << nom3 << endl
	   << nom4 << endl ;
    }

 

  string nom ;

  if ( FileExists(nom1) )
    {
      nom = nom1 ;
    }
  else  
    {
      if ( FileExists(nom2) )
	{
	  nom = nom2 ;
	}
      else  
	{
	  if ( FileExists(nom3) )
	    {
	      nom = nom3 ;
	    }
	  else  
	    {
	      if ( FileExists(nom4) )
		{
		  nom = nom4 ;
		}
	      else
		{
		  cerr << "No existing file " << endl ; 
		  exit(0);
		}
	    }
	}
    }
    
  FitsHeader header(nom.c_str());
  int taille_x = header.KeyVal("NAXIS1") ;
  int taille_y = header.KeyVal("NAXIS2") ;

  cerr << "Taille image entree " << taille_x << " " << taille_y << endl ;

  
  if ( !FileExists(nom1) )
    {
      cerr << "Can't find " << nom1 << endl ;
      //exit(0);
    }

  if ( !FileExists(nom2)  )
    {
      cerr << "Can't find " << nom2 << endl ;
      //exit(0);
    }


  if ( !FileExists(nom3)  )
    {
      cerr << "Can't find " << nom3 << endl ;
      //exit(0);
    }


  if ( !FileExists(nom4)  )
    { 
      cerr << "Can't find " << nom4 << endl ;
      //exit(0);
    }









  int taillex_scale = Taille_Scale(taille_x, delta);
  int tailley_scale = Taille_Scale(taille_y, delta);


  int Tx = tailley_scale + taillex_scale ;
  int Ty = 3 * taillex_scale ;
  Image mosa(Tx, Ty);

  if ( FileExists(nom1)  )
      {
       FitsImage *pim1 = new FitsImage(nom1.c_str());      
       cerr << "Tx, Ty : " << pim1->Nx() << " " << pim1->Ny() << endl ;
       
       Image * pim1s = Rescale(*pim1, delta);
       cerr << "Txs, Tys : " << pim1s->Nx() << " " << pim1s->Ny() << endl ;
       delete pim1;
 
 
       Image *pim1r = Retourne(pim1s);
       delete pim1s ;
       cerr << "Tx, Ty : " << pim1r->Nx() << " " << pim1r->Ny() << endl ;
  
       cerr << "Adding ccd 1 " << endl ;
       Add(*pim1r,  mosa, 0, 0) ;
       delete pim1r ;
      }




   if ( FileExists(nom2)  )
      {
	  FitsImage *pim2 = new FitsImage(nom2.c_str());  
	  Image * pim2s = Rescale(*pim2, delta);
	  delete pim2 ;
	  Image *pim2r = Retourne5(pim2s);
	  delete pim2s ;
	  cerr << "Adding ccd 2 " << endl ;
	  Add(*pim2r,  mosa,  tailley_scale, 0) ;
	  delete pim2r ; 
      }


   if ( FileExists(nom3)  )
      {
       FitsImage *pim3 = new FitsImage(nom3.c_str());  
       Image * pim3s = Rescale(*pim3, delta);
       delete pim3 ;
       Image *pim3r = Retourne(pim3s);
       delete pim3s ;
       cerr << "Adding ccd 3 " << endl ;
       Add(*pim3r,  mosa,  0, 2*taillex_scale) ;
       delete pim3r ;
      }


   if ( FileExists(nom4)  )
      {
       FitsImage *pim4 = new FitsImage(nom4.c_str());  
       Image * pim4s = Rescale(*pim4, delta);
       delete pim4 ;
       Image *pim4r = Retourne(pim4s);
       delete pim4s ;
       cerr << "Adding ccd 4 " << endl ;
       Add(*pim4r,  mosa, 0,  taillex_scale) ;
       delete pim4r ;
      }

 /* on a maintenant des images:

        ___________
        |         |
        |   3     |
        o----------
        ----------- ------
        |         | |    |
        |   4     | |    |
        o---------- |5(2)|
        ----------- |    |
taille_x|         | |    |
        |   1     | |    |
        o---------- o-----
          taille_y
*/
 

  
  FitsImage("mosa.fits", mosa);

}

int
Taille_Scale(int Taille, int delta)
{ 
    float scale = 1. / delta ;
    int T2 = (int) ( Taille*scale - 1.) ;
    return(T2);
}


Image *
Rescale(Image & image, int delta)
{
  float scale = 1. / delta ;
  float scale2 = scale * scale ;
  int Tx = Taille_Scale(image.Nx(), delta) ;
  int Ty = Taille_Scale(image.Ny(), delta) ;

  Image * pim = new Image(Tx,Ty) ;
  for (int i = 0 ; i < Tx ; i++)    
    for (int j = 0 ; j < Ty ; j++)
      {
	double u = 0. ;
	for (int k = 0; k < delta ; k++)	  
	  for (int l = 0; l < delta ; l++)
	    u += image(delta * i + k , delta * j  + l ) ;
	u *= scale2 ;
	(*pim)(i,j) = u ;
      }
  return(pim);
}


void
Add(Image const & image1, Image & imagetot, int xor_tot, int yor_tot)
{
 for (int i = 0 ; i <  image1.Nx() ; i++)    
   for (int j = 0 ; j < image1.Ny() ; j++)
     {
      imagetot(  xor_tot + i ,  yor_tot + j ) = image1(i,j);
     }
}



Image *
Retourne(Image *pim1)
{
  Image * pim1r = new Image(pim1->Ny(), pim1->Nx());  
  cerr << "Tx, Ty : " << pim1r->Nx() << " " << pim1r->Ny() << endl ;

  for (int i = 0 ; i <  pim1->Nx() ; i++)    
   for (int j = 0 ; j <  pim1->Ny() ; j++)
     {
       int yy = pim1->Nx() -1 -i ;
       int xx = j ;
					    
       (*pim1r)(xx,yy) = (*pim1)(i,j);
     }
 return(pim1r);

}


Image *
Retourne5(Image *pim1)
{
  Image * pim1r = new Image(pim1->Nx(), pim1->Ny());
  cerr << "Tx, Ty : " << pim1r->Nx() << " " << pim1r->Ny() << endl ;
  for (int i = 0 ; i <  pim1->Nx() ; i++)    
   for (int j = 0 ; j <  pim1->Ny() ; j++)
     {
       int xx = pim1->Nx() -1 -i ;
       int yy = pim1->Ny() -1 -j ;
					    
       (*pim1r)(xx,yy) = (*pim1)(i,j);
     }
 return(pim1r);
}

