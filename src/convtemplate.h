// This may look like C code, but it is really -*- C++ -*-



static void 
Convole_LDEMICOTE(Image const &  imgref , 
Image const & noyau , Image & imgconv )
/*convolution de imgref, resultat ds imgconv. on convole avec le  
noyau non normalise . ilpeut y a avoir un masque (=NULL) par defaut.
Ds ce cas, on copie imgref ds imgconv et on ne convole que la ou il y a 1
su le masque */
{
      TRY{
    int i ,  j , l , k , d = LDEMICOTE;  ;
    int tx , ty ;
    const Pixel  *penhautagauche;
    Pixel *parrivee;
    int deltax ;
    double Norma = 0 ;
    for(i=0 ; i<2*d + 1 ; i++)
      for (j=0 ; j<2*d + 1 ; j++)
	{ 
	  Norma += noyau(i,j)   ;
	}


    tx = imgref.Nx();
    ty = imgref.Ny();

    penhautagauche = imgref.begin() ;
    parrivee = imgconv.begin() ;

    parrivee += d*tx + d ;

    deltax = tx - 2* d - 1;

    for(i=0; i < ty - 2* d ; i++ )
      { 
        for(j=0; j < tx - 2* d ; j++ )
	  { 
	    register const Pixel  * pimage = penhautagauche ;
	    register double S=0;
	    register const Pixel  * pnoyau =noyau.begin() ;
	    for (k= 0; k < 2* d +1  ; k++ )
	      { 
		for (l= 0; l < 2*d + 1 ; l++ )  
		  { 
		    Pixel f = *pimage ;
		    double  nn = *pnoyau  ;
		    S += (double) ( nn * f ) ;
		    pimage++ ;
		    pnoyau++ ;
		    //S += *pimage++ * *pnoyau++ ;
		  }
		pimage += deltax;
	      } 
	    double r = S / Norma ;
	    *parrivee++ = (Pixel) r ;
	    penhautagauche++ ;
	  }
	penhautagauche += 2* d   ; 
	parrivee += 2*d ;
      } 
    }
    CATCHALL{THROW_SAME;}
    ENDTRY
    return ;
}



static void 
Convole_LDEMICOTE(Image const & imgref ,   Image const & noyau,
		  Image const  & mask , Image & imgconv)
/*convolution de imgref, resultat ds imgconv. on convole avec le 
noyau non normalise 
Ds ce cas, on copie imgref ds imgconv et on ne convole que la ou il y a 1
su le masque */
{
  //TIMEF
   TRY{
    int i ,  j , l , k , d = LDEMICOTE;  ;
    int tx , ty ;
    const Pixel  *penhautagauche;
     Pixel  *parrivee;
    int deltax ;

    double Norma = 0 ;
    for(i=0 ; i<2*d + 1 ; i++)
      for (j=0 ; j<2*d + 1 ; j++)
	{ 
	  Norma += noyau(i,j)   ;
	}



    tx = imgref.Nx();
    ty = imgref.Ny();

    penhautagauche = imgref.begin() ;
    parrivee = imgconv.begin() ;

    parrivee += d*tx + d ;

    deltax = tx - 2* d - 1;




    /* ne vont etre convoles que certains pixels , a 1 sur le masque.
       les autres restent tels quels */
    imgconv =  imgref ; 
    const Pixel* pmask ;
    pmask= mask.begin();
    pmask += d*tx + d ;
    /* mis sur le premier pixel qui peut etre convole */
    // le mask est en avance sur imgref et imgconv

    for(i=0; i < ty - 2* d ; i++ )
      { 
        for(j=0; j < tx - 2* d ; j++ )
	  { 
	    register const Pixel * pimage = penhautagauche ;
	    register double S=0; 
	    register const Pixel * pnoyau = noyau.begin() ;
	    if ( *pmask  ) /* convolution */
	      {
		for (k= 0; k < 2* d +1  ; k++ )
		  { 
		    for (l= 0; l < 2*d + 1 ; l++ )  
		      { 
			S += *pimage++ * *pnoyau++ ;
		      }
		    pimage += deltax;
		  } 
		*parrivee++ =  (Pixel) (S / Norma );
	      }
	    else
	      parrivee++ ; // on ne fait rien
           
	    pmask++ ;
	    penhautagauche++ ;
	  }
        penhautagauche += 2* d   ; 
        parrivee += 2 * d ;
        pmask += 2 * d ;
      } 
    }
    CATCHALL{THROW_SAME;}
    ENDTRY
    return; 
}

static void Convole_LDEMICOTE(Image const & imgref ,   double s ,
  Image const & mask ,   Image & imgconv)
{
  TRY{
    Image noyau = Convole_MakeNoyau(s,LDEMICOTE);   
    Convole_LDEMICOTE(imgref, noyau , mask , imgconv);
    } 
CATCHALL{THROW_SAME;}
ENDTRY
}

void Convole_LDEMICOTE(Image const & imgref ,   double s ,
    Image & imgconv)
{
  TRY{
  Image noyau = Convole_MakeNoyau(s,LDEMICOTE);   
  Convole_LDEMICOTE(imgref, noyau , imgconv);
  } 
  CATCHALL{THROW_SAME;}
ENDTRY
}



static void Convole_LDEMICOTE(Image const & imgref ,   double sx , double sy ,
    Image & imgconv)
{
  TRY{
  Image noyau = Convole_MakeNoyau(sx , sy,LDEMICOTE);   
  Convole_LDEMICOTE(imgref, noyau , imgconv);
} 
  CATCHALL{THROW_SAME;}
  ENDTRY
}


static void Convole_LDEMICOTE(Image const & imgref ,   double sx , double sy ,
   Image const  & mask , Image & imgconv)
{
  TRY{
  Image noyau = Convole_MakeNoyau(sx, sy,LDEMICOTE);   
  Convole_LDEMICOTE(imgref, noyau , mask ,imgconv);
  } 
  CATCHALL{THROW_SAME;}
  ENDTRY
}




/**********************************************/
/******** CONVOLUTION 1D **********************/
/**********************************************/

static void 
Convole1DX_LDEMICOTE(Image const &  imgref , double * noyau , 
		     Image & imgconv )
{
  //convolution de imgref, resultat ds imgconv. on convole avec le  
  //noyau non normalise .  

  //TIMEF 
    TRY{
    int i ,  j ,  k , d = LDEMICOTE;  ;
    int tx , ty ;
    const Pixel  *penhautagauche;
    Pixel *parrivee;
    double Norma = 0 ;
    for(i=0 ; i<2*d + 1 ; i++)
      { 
	Norma += noyau[i]  ;
      }

    tx = imgref.Nx();
    ty = imgref.Ny();

    penhautagauche = imgref.begin() ;
    parrivee = imgconv.begin() ;

    parrivee +=  d ;
    int Ntot = 0 ;
    for(i=0; i < ty  ; i++ )
      { 
        for(j=0; j < tx - 2* d ; j++ )
	  { 
	    register const Pixel  * pimage = penhautagauche ;
	    register double S=0; 
	    register double * pnoyau = noyau ;
	    for (k= 0; k < 2* d +1  ; k++ )
	      {
		double nn = *pnoyau ;
		Pixel f  = *pimage ;
		double ff = (double) f ;
		S += nn * ff ;
		pimage++ ;
		pnoyau++ ; 
	      } 
	    Ntot++ ;
	    *parrivee++ = (Pixel) (S / Norma ) ;
	    penhautagauche++ ;
	  }
	penhautagauche += 2* d   ; 
	parrivee += 2*d ;
      }
    } 
    CATCHALL{THROW_SAME;}
    ENDTRY
    return ;
}

static void 
Convole1DX_LDEMICOTE(Image const & imgref , double * noyau , 
		     Image const  & mask , Image & imgconv)
{
  //convolution de imgref, resultat ds imgconv. on convole avec le 
  //noyau non normalise 
  // Ds ce cas, on copie imgref ds imgconv et on ne convole que la ou il y a 1
  // sur le masque 

  //TIMEF
  TRY{
    int i ,  j ,  k , d = LDEMICOTE;  ;
    int tx , ty ;
    const Pixel  *penhautagauche;
    Pixel *parrivee;

    double Norma = 0 ;
    for(i=0 ; i<2*d + 1 ; i++)
      { 
	Norma += noyau[i]  ;
      }




    tx = imgref.Nx();
    ty = imgref.Ny();

    penhautagauche = imgref.begin() ;
    parrivee = imgconv.begin() ;

    parrivee +=  d ;




    /* ne vont etre convoles que certains pixels , a 1 sur le masque.
       les autres restent tels quels */
    imgconv = imgref ; 
    const Pixel * pmask ;
    pmask= mask.begin();
    pmask += d ; /* mis sur le premier pixel qui peut etre convole */

    for(i=0; i < ty  ; i++ )
      { 
        for(j=0; j < tx - 2* d ; j++ )
	  { 
	    register const Pixel * pimage = penhautagauche ;
	    register double S=0; 
	    register double * pnoyau = noyau ;
	    if ( *pmask  ) /* convolution */
	      {
		for (k= 0; k < 2* d +1  ; k++ )
		  { 
		    S += *pimage++ * *pnoyau++  ;   
		  } 
		*parrivee++ = (Pixel) (S / Norma ) ;
                      

	      }
	    else
	      parrivee++ ; // on ne fait rien
           
	    pmask++ ;
	    penhautagauche++ ;
	  }
        penhautagauche += 2* d   ; 
        parrivee += 2 * d ;
        pmask += 2 * d ;
      }
    }  
    CATCHALL{THROW_SAME;}
    ENDTRY
  
    return; 
}

static void 
Convole1DY_LDEMICOTE(Image const &  imgref , 
		     double * noyau , Image & imgconv )
  // convolution de imgref, resultat ds imgconv. on convole avec le  
  // noyau non normalise .
{
  ////TIMEF
  TRY{
    int i ,  j ,  k , d = LDEMICOTE;  ;
    int tx , ty ;
    const Pixel  *penhautagauche;
    Pixel *parrivee;
    double Norma = 0 ;
    for(i=0 ; i<2*d + 1 ; i++)
      { 
	Norma += noyau[i]  ;
      }


    tx = imgref.Nx();
    ty = imgref.Ny();

    penhautagauche = imgref.begin() ;
    parrivee = imgconv.begin() ;


    parrivee +=  d  * tx;
    for(i=0; i < ty - 2* d  ; i++ )
      { 
        for(j=0; j < tx  ; j++ )
	  { 
	    register const Pixel  * pimage = penhautagauche ;
	    register double S=0; 
	    register double * pnoyau = noyau ;
	    for (k= 0; k < 2* d +1  ; k++ )
	      {
		Pixel f = *pimage ;
		double ff = (double) f ;
		//double nn = *pnoyau ;
		double nn = noyau[k] ;
		S += ff * nn ;   
		pimage +=  tx ; pnoyau++  ;
	      }
	    double r = S / Norma ;
	    *parrivee++ = (Pixel) r ;
	    penhautagauche++ ;
	  }
       
      } 
    }  
    CATCHALL{THROW_SAME;}
    ENDTRY
    return ;
}

static void 
Convole1DY_LDEMICOTE(Image const & imgref ,   double * noyau,
		     Image const  & mask , Image & imgconv)
{
  //  convolution de imgref, resultat ds imgconv. on convole avec le 
  //noyau non normalise 
  //Ds ce cas, on copie imgref ds imgconv et on ne convole que la ou il y a 1
  //sur le masque 
  //TIMEF
    TRY{
    int i ,  j ,  k , d = LDEMICOTE;  ;
    int tx , ty ;
    const Pixel  *penhautagauche;
    Pixel *parrivee;



    double Norma = 0 ;
    for(i=0 ; i<2*d + 1 ; i++)
      { 
	Norma += noyau[i]  ;
      }



    tx = imgref.Nx();
    ty = imgref.Ny();

    penhautagauche = imgref.begin() ;
    parrivee = imgconv.begin() ;

    parrivee += d*tx  ;

    // ne vont etre convoles que certains pixels , a 1 sur le masque.
    //les autres restent tels quels 
    imgconv =  imgref ; 
    const Pixel * pmask ;
    pmask= mask.begin();
    pmask += d*tx  ; /* mis sur le premier pixel qui peut etre convole */

    for(i=0; i < ty - 2* d ; i++ )
      { 
        for(j=0; j < tx  ; j++ )
	  { 
	    register const Pixel * pimage = penhautagauche ;
	    register double S=0; 
	    register double * pnoyau = noyau ;
	    if ( *pmask ) /* convolution */
	      {
		for (k= 0; k < 2* d +1  ; k++ )
		  {
		    S += *pimage * *pnoyau  ;   
		    pimage +=  tx ; pnoyau++  ;
		  } 
		*parrivee++ = (Pixel) (S / Norma ) ;
                    
	      }
	    else
	      parrivee++ ; // on ne fait rien
           
	    pmask++ ;
	    penhautagauche++ ;
	  }
      } 
 
    }  
    CATCHALL{THROW_SAME;}
    ENDTRY 
    return; 
}

static void 
Convole1DX_LDEMICOTE(Image const & imgref ,   double s ,
		     Image const & mask ,   Image & imgconv)
{
  double * noyau;
  TRY{
  noyau = Convole_MakeNoyau_1D(s,LDEMICOTE);   
  Convole1DX_LDEMICOTE(imgref, noyau , mask , imgconv);
  }
   CATCHALL{THROW_SAME;}
  ENDTRY 
    return; 
}

static void 
Convole1DX_LDEMICOTE(Image const & imgref ,   double s ,
		     Image & imgconv)
{
  double * noyau;
  TRY{
  noyau = Convole_MakeNoyau_1D(s,LDEMICOTE);    
  Convole1DX_LDEMICOTE(imgref, noyau , imgconv);
  }
  CATCHALL{THROW_SAME;}
  ENDTRY 
    return; 
}

static void 
Convole1DY_LDEMICOTE(Image const & imgref ,   double s ,
		     Image const & mask ,   Image & imgconv)
{
  double * noyau;
  TRY{
  noyau = Convole_MakeNoyau_1D(s,LDEMICOTE);   
  Convole1DY_LDEMICOTE(imgref, noyau , mask , imgconv);
  }
  CATCHALL{THROW_SAME;}
  ENDTRY 
    return; 
}

static void 
Convole1DY_LDEMICOTE(Image const & imgref ,   double s ,
		     Image & imgconv)
{
  double * noyau;
  TRY{
  noyau = Convole_MakeNoyau_1D(s,LDEMICOTE);    
  Convole1DY_LDEMICOTE(imgref, noyau , imgconv);
  }
  CATCHALL{THROW_SAME;}
    ENDTRY 
    return; 
}

