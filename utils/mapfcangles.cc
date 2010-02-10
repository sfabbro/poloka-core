// -*- C++ -*-
// 
#include <math.h>
#include <iostream>

#include <frame.h>
#include <gtransfo.h>

#define CHECK_BOUNDS
#include <fitsimage.h>
#include <imageutils.h>
#include <wcsutils.h>




#include <list>

struct AstromData
{
  CountedRef<Gtransfo> wcs;
  CountedRef<Gtransfo> inv_wcs;
  Frame radecFrame;
  Point center;
  int chip;

  AstromData(const FitsHeader &Head)
  {
    chip = Head.KeyVal("TOADCHIP");
    Gtransfo *wcs_temp = NULL;
    WCSFromHeader(Head, wcs_temp);
    wcs = wcs_temp;
    //radecFrame = pixFrame for now
    radecFrame=Frame(Head);
    inv_wcs = wcs->InverseTransfo(0.1,radecFrame); 
    // becomes what it should be;
    radecFrame = ApplyTransfo(radecFrame, *wcs);


  }

  const Gtransfo *Wcs() const { return wcs;}
  const Gtransfo *InvWcs() const { return inv_wcs;}
    
};

typedef list<AstromData> adList;



// in the ccd order : CCD0 -> CCD35)
static int xstart[4][18] = 
  {
    {132, 197, 264 ,329, 396, 461, 528, 593, 661, 726, 793, 858, 925, 990, 1058, 1123, 1190, 1255},
    {132, 197, 264, 329, 396, 461, 529, 593, 661, 726, 793, 858, 926, 990, 1058, 1123, 1190, 1255},
    {132, 197, 264, 329, 396, 461, 529, 594, 661, 726, 793, 858, 925, 990, 1058, 1123, 1190, 1255},
    {132, 197, 264, 329, 396, 461, 529, 594, 661, 726, 793, 858, 926, 991, 1058, 1123, 1190, 1255}
  };

static int ybottom[4][18] = 
  {{921, 921, 921, 921, 921, 921, 921, 921, 921, 921, 921, 921, 921, 921, 921, 921, 921, 921},
   {606,606 ,606 ,606 ,606 ,606 ,606 ,606 ,607 ,607 ,606 ,606 ,606 ,606 ,607 ,607 ,606 ,606},
   {315, 315, 315, 315, 315, 315, 315, 315, 315, 315, 315, 315, 315, 315, 315, 315, 315, 315},
   {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0}
  };


//int resamp_factor = 16;

static const int nx_prev = 1452;
static const int ny_prev=1210;


/* routine qui retourne un numero d'ampli (0:71) a partir de coordonnees dans le "preview".
   le numero de la CCD, c'est WhichHalfCCD()/2.
   Attention: selon la moitie de la camera consideree, l'ampli correspondant a x<1024 est le premier
   ou le 2e des 2 amplis qui constituent une CCD. 
   Les CCDs<18 sont retournes ((x,y)=(0,0) est en haut a droite).
*/
int WhichHalfCCD(const int xPrev, const int yPrev)
{
  if (xPrev<0 || xPrev>= nx_prev || yPrev<0 || yPrev >= ny_prev) return 99;
  int xSize = 64;
  int ySize = 289;
  int xOffset = xstart[0][0];
  int i0_guess = (xPrev-xOffset)/xSize;
  int j0_guess = 3-(yPrev/ySize);
  static const int offsets[3] = {0,-1,1};
  for (int ik=0; ik<3; ++ik)
    {
      int i0=i0_guess+offsets[ik];
      if (i0<0 || i0>=18) continue;
      for (int jk=0;jk<3; ++jk)
	{
	  int j0 = j0_guess+offsets[jk];
	  if (j0<0 || j0>=4) continue;
	  if ((xPrev >= xstart[j0][i0] && xPrev<xstart[j0][i0]+xSize)
	      && (yPrev >= ybottom[j0][i0] && yPrev<ybottom[j0][i0]+ySize)
	      )
	    return 18*j0+i0;
	}
    }
  return 99;
}




int main(int nargs, char **args)
{
  FitsImage ax("ax.fits",nx_prev,ny_prev);
  FitsImage ay("ay.fits",nx_prev,ny_prev);
  FitsImage theta("theta.fits",nx_prev,ny_prev);
  FitsImage phi("phi.fits",nx_prev,ny_prev);

  ax.ModKey("BITPIX",-32);
  ay.ModKey("BITPIX",-32);
  theta.ModKey("BITPIX",-32);
  phi.ModKey("BITPIX",-32);
    
  FitsHeader toto(args[1]);
  Point center(double(toto.KeyVal("RA_DEG")), double(toto.KeyVal("DEC_DEG")));  
  TanRaDec2Pix radec2pix( GtransfoLin(), center );  
  
  for (int i=1; i< nargs; ++i)
    {
      const char *arg = args[i];
      FitsHeader head(arg);
      Gtransfo *wcs_temp = NULL;
      WCSFromHeader(head, wcs_temp);
      const Gtransfo *wcs = wcs_temp;
      int ccd = head.KeyVal("TOADCHIP");
      double cd1_1 = head.KeyVal("CD1_1");
      int reversed_x = (cd1_1>0) ? 1 : 0;
      double cd2_2 = head.KeyVal("CD2_2");
      int reversed_y = (cd2_2>0) ? 0 : 1;

      int j0 = ccd/9;
      int col0 = ccd-9*j0;

      for (int i=8; i<2048; i+=16)
	{
	  int iamp;
	  if (reversed_x)
	    iamp = (i>1024) ? 0 : 1;
	  else
	    iamp = (i>1024) ? 1 : 0;
	  int i_on_amp = i%1024;
	  int xprev;
	  if (reversed_x)
	    xprev = xstart[j0][2*col0+iamp]+63-(i_on_amp-8)/16;
	  else
	    xprev = xstart[j0][2*col0+iamp]+(i_on_amp-8)/16;
	  for (int j=8; j< 4612; j+=16)
	    {
	      int yprev;
	      if (reversed_y)
		yprev = ybottom[j0][2*col0+iamp]+288-(j-8)/16;
	      else
		yprev = ybottom[j0][2*col0+iamp]+(j-8)/16;
	      Point where_pix(i, j);
	      Point where_sky = wcs->apply(where_pix);
	      Point where_tp;
	      radec2pix.apply( where_sky.x, where_sky.y, 
			       where_tp.x,  where_tp.y );
	      
	      ax(xprev,yprev) = where_tp.x; // -center.x;
	      ay(xprev,yprev) = where_tp.y; // center.y;
	      theta( xprev, yprev ) = sqrt( where_tp.x*where_tp.x + where_tp.y*where_tp.y ) * M_PI / 180.;
	      phi( xprev, yprev )   = -atan2( where_tp.y, where_tp.x );

	      /*
		ax(xprev,yprev) = 2*col0+iamp;
		ay(xprev,yprev) = j0;
	      */
	    }
	}
    } // end loop on input files
  return EXIT_SUCCESS;
}
