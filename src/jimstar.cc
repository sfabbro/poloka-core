#include <iostream>
#include <fstream>
#include <math.h>
#include <cmath>


#include "sestar.h"
#include "usnoutils.h"
#include "wcsutils.h"
#include "frame.h"
#include "starmatch.h"
#include "listmatch.h"
#include "jimstar.h"
#include "astroutils.h"

#ifndef M_PI
#define  M_PI           3.14159265358979323846  /* pi */
#endif



JimStar::JimStar()
  : BaseStar(0.,0.,0.) {
  Set_to_Zero();
}

JimStar::JimStar(double xx, double yy, double ff)
  : BaseStar(xx,yy,ff) {
  Set_to_Zero();
}

void JimStar::Set_to_Zero() {
  name="";
  ra=0;
  dec=0;
  imag=0.0;
  rmag=0.0;
  gmag=0.0;
  zmag=0.0;
  cgal=0.0;
}

void JimStar::Read(istream& r, const char *Format) {
  //  char Name[16];
  //  char Ra[12];
  //  char Dec[12];
  int format = DecodeFormat(Format, "JimStar");
  
  r >> x;
  r >> y;
  r >> imag;
  r >> rmag;
  r >> gmag;
  r >> zmag;
  r >> cgal;

  ra = x;
  dec = y;
}


JimStar* JimStar::read(istream& r, const char *Format) {
  
  JimStar *pstar = new JimStar();
  pstar->Read(r,Format);
  return(pstar);
}

void JimStar::dumpn(ostream& s) const {
  s << " x : " << x;
  s << " y : " << y;
  s << " ra : " << ra;
  s << " dec : " << dec;
  s << " imag : " << imag;
  s << " rmag : " << rmag;
  s << " gmag : " << gmag;
  s << " zmag : " << zmag;
  s << " cgal : " << cgal;
}

void JimStar::dump(ostream& s) const {
  dumpn(s);
  s << endl;
}

void JimStar::writen(ostream& s) const {
  s << x << " " ;
  s << y << " " ;
  s << ra << " " ;
  s << dec << " " ;
  s << imag << " " ;
  s << rmag << " " ;
  s << gmag << " " ;
  s << zmag << " " ;
  s << cgal << " " ;
}

void JimStar::write(ostream& s) const {
  writen(s);
  s << endl; 
}

string JimStar::WriteHeader_(ostream & pr, const char *i) const {
  
  if (i==NULL) i="";
  //string baseStarFormat =  BaseStar::WriteHeader_(pr, i);
  pr 
    << "# x"<< i << " : x...  " << endl
    << "# y"<< i << " : y... " << endl
    << "# ra"<< i << " : ra...  " << endl
    << "# dec"<< i << " : dec... " << endl
    << "# imag"<< i << " : std mag in I " << endl
    << "# rmag"<< i << " : std mag in R " << endl
    << "# gmag"<< i << " : std mag in G " << endl
    << "# zmag"<< i << " : std mag in Z " << endl
    << "# cgal"<< i << " : cgal " << endl;
    

  //static char format[256];
  string format = " JimStar 1";
  //sprintf(format,"%s JimStar %d",baseStarFormat, 1);
  return format;
}

JimStarList* GetSelectedJimStarList(const FitsHeader &header, const Frame &W, const string &standardfile)
{
  Gtransfo* Pix2RaDec=0;
  WCSFromHeader(header, Pix2RaDec);
  Gtransfo *RaDec2Pix = Pix2RaDec->InverseTransfo(0.01,W);
  
  Frame radecW = (W.ApplyTransfo(*Pix2RaDec)).Rescale(1.1);
  // cout << "radecW = " << radecW << endl;
  
  JimStarList *stdstarlist = new JimStarList(standardfile);
  
  if (!stdstarlist || stdstarlist->size()==0)
    {
      cout << "Bad standard file !!!" << endl;
      return NULL;
    }
  // test
  
  // return stdstarlist;

  for (JimStarIterator si = stdstarlist->begin(); si != stdstarlist->end(); ) {
    JimStar &pstar = **si;
    if(!radecW.InFrame(pstar)) {
      si = stdstarlist->erase(si);
      continue;
    }
    //cout << "b = " << pstar << endl;
    RaDec2Pix->apply(pstar,pstar);
    //cout << "a = " << pstar << endl;
    if (W.InFrame(pstar))
      ++ si;
    else
      si = stdstarlist->erase(si);
  }
  return stdstarlist;
}


/************************** FINDEFINITION JimStar ************************/

#include "starlist.h"
#include "starlist.cc" /* since starlist is a template class */

template class StarList<JimStar>;

