// This may look like C code, but it is really -*- C++ -*-

#ifndef DAOPHOTIO__H
#define DAOPHOTIO__H

#include <iomanip>
#include <fstream>
#include "daophotpsf.h"
#include "sestar.h"
#include "reducedimage.h"

//! the zeropoint used by DAOPHOT when computing instrumental magnitudes
const double DAOPHOT_APER_ZP = 25.;

//! the shift in pixels between DAOPHOT coordinate system and TOADS one.
const double DAOPHOT_TOADS_SHIFT = -1.;

//! the type of star list in the daophot private style files
int DaoFileNumber(const DbImageCatalogKind FileType);

//! reads a SEStar from a stream, according to the filetype
template<DbImageCatalogKind filetype> 
void read_dao_star(istream &daostream, SEStar &star);

//! reads and fills in the star list from a stream given a filetype
template<DbImageCatalogKind filetype> 
void read_dao_starlist(istream &daostream, SEStarList &Stars, const double zp=DAOPHOT_APER_ZP)
{
  char c;
  while (daostream >> c) // end-of-file test
    {
      daostream.unget();
      SEStar *star = new SEStar;
      double x,y, mag;
      daostream >> star->N() >> x >> y >> mag;
      // DAOPHOT to TOADS conversion x,y and flux
      star->x = x + DAOPHOT_TOADS_SHIFT; 
      star->y = y + DAOPHOT_TOADS_SHIFT;
      star->flux = pow(10.,(0.4*(zp - mag)));
      read_dao_star<filetype>(daostream, *star); 
      Stars.push_back(star);
    }
}

//! reads and fills in the star list and the header from a stream given a filetype and the name of the file
template<DbImageCatalogKind filetype> 
void read_dao(const string &FileName, SEStarList &Stars, const double zp=DAOPHOT_APER_ZP)
{
  ifstream daostream(FileName.c_str());
  if (!daostream) 
    {
      cerr << " read_dao : unable to open file '" 
	   << FileName << "'" << endl;
      return;
    }
  string dummy;
  getline(daostream, dummy);
  getline(daostream, dummy);
  read_dao_starlist<filetype>(daostream, Stars, zp); 
}

//! reads and fills in the star list from a stream given a filetype and a DbImage
template<DbImageCatalogKind filetype> 
void read_dao(const DbImage& Dbim, SEStarList &Stars)
{
  double zp = DAOPHOT_APER_ZP;
  if (((filetype == DaophotAls)||
       (filetype == DaophotPk) ||
       (filetype == DaophotNst))
      && FileExists(Dbim.ImagePsfName()))
    {
      DaoPsf psf(Dbim);
      zp = psf.PsfMag();
    }
  read_dao<filetype>(Dbim.ImageCatalogName(filetype), Stars, zp);
}

//! writes a SEStar to a stream, according to the filetype
template<DbImageCatalogKind filetype> 
void write_dao_star(ostream &daostream, const SEStar &star);

//! writes a properly formatted daophot style header
template<DbImageCatalogKind filetype> 
void write_dao_header(ostream &daostream, const int Nx, const int Ny, 
		      const double LowBad, const double HighBad, 
		      const double Threshold, const double Aperture, 
		      const double Gain, const double ReadNoise, const double FitRad)
{
  int nl = DaoFileNumber(filetype);
  daostream << " NL   NX   NY  LOWBAD HIGHBAD  THRESH     AP1  PH/ADU  RNOISE    FRAD\n"; 
  int ndigits=0;
  if (HighBad > 99999.9) ndigits=0;
  daostream << setiosflags(ios::fixed);
  daostream << " " 
	    << setw(2) << nl 
	    << setw(5) << Nx 
	    << setw(5) << Ny 
	    << setw(8) << setprecision(1) << LowBad 
	    << setw(8) << setprecision(ndigits) << HighBad 
	    << setw(8) << setprecision(2) << Threshold 
	    << setw(8) << Aperture << setw(8) << Gain << setw(8) << ReadNoise << setw(8) << FitRad;
  daostream << endl << endl;
  daostream << resetiosflags(ios::fixed);
}

//! writes a daophot header from option style daophot to be consistent
template<DbImageCatalogKind filetype> 
void write_dao_header(ostream &daostream, const int ncol, const int nrow, 
		      const float sky, const float sigma, const float *opt)
{
  write_dao_header<filetype>(daostream, ncol, nrow, sky-opt[2]*sigma, 
			     opt[3], opt[5]*sigma, opt[4], opt[1], opt[0], opt[11]);
}

//! writes a properly formatted daophot style header from a ReducedImage
template<DbImageCatalogKind filetype> 
void write_dao_header(ostream &daostream, const ReducedImage &Rim)
{
  write_dao_header<filetype>(daostream, Rim.XSize(), Rim.YSize(),
			     Rim.BackLevel() - 7* Rim.SigmaBack(),
			     Rim.Saturation(), 4*Rim.SigmaBack(), 
			     Rim.Seeing()*2.3548, Rim.Gain(), Rim.ReadoutNoise(), 
			     Rim.Seeing()*2.3548);
}


//! writes a SEStarList into a stream given a Daophot filetype
template<DbImageCatalogKind filetype> 
void write_dao_starlist(ostream &daostream, const SEStarList &Stars, const double zp =DAOPHOT_APER_ZP)
{
  daostream << setiosflags(ios::fixed);
  for (SEStarCIterator it=Stars.begin(); it!=Stars.end(); ++it)
    {
      const SEStar *star = *it ;
      double mag = (star->flux > 0)? -2.5*log10(star->flux) + zp : zp;
      daostream << setw(6) << star->N()
		<< setprecision(3) 
	// DAOPHOT-TOADS conventions
		<< setw(9) << star->x - DAOPHOT_TOADS_SHIFT 
		<< setw(9) << star->y - DAOPHOT_TOADS_SHIFT
		<< setw(9) << mag;
      write_dao_star<filetype>(daostream, *star);
      daostream << endl;
    }
  daostream << resetiosflags(ios::fixed);
}

//! writes a SEStarList for a ReducedImage with proper file type
template<DbImageCatalogKind filetype> 
void write_dao(const ReducedImage &Rim, const SEStarList &Stars)
{
  ofstream daostream((Rim.ImageCatalogName(filetype)).c_str());
  write_dao_header<filetype>(daostream, Rim);
  write_dao_starlist<filetype>(daostream, Stars); 
}

//! writes a SEStarList in a filename. Needs some more info from the image.
template<DbImageCatalogKind filetype> 
void write_dao(const string FileName, const int ncol, const int nrow, 
	       const float sky, const float sigma, 
	       const float *opt, const SEStarList &Stars)
{
  ofstream daostream(FileName.c_str());
  write_dao_header<filetype>(daostream, ncol, nrow, sky, sigma, opt);
  write_dao_starlist<filetype>(daostream, Stars); 
}


//! combines a DAOPHOT and SExtractor file of a DbImage into a SEStarList
template<DbImageCatalogKind filetype> 
void CombineSextractorDaophot(const DbImage &DbIm, 
			      const bool ReplacePos=true, 
			      const bool ReplaceFlux=false,
			      const bool ReplaceSky=false)
{
  if (filetype == SExtractor) 
    {
      cout << " CombineSextractorDaophot: WARNING : combining the same sextractor catalogs !!!" << endl;
      return;
    }

  SEStarList sexstars(DbIm.ImageCatalogName(SExtractor));  
  SEStarList daostars;
  read_dao<filetype>(DbIm, daostars);

  SEStarIterator enddao = daostars.end();
  size_t ndeleted = 0;

  for (SEStarIterator itsex = sexstars.begin(); itsex != sexstars.end(); ++itsex)
    {
      SEStarCIterator itdao = daostars.begin();
      while (itdao != enddao && (*itdao)->N() != (*itsex)->N()) ++itdao;
      if (itdao == enddao) 
	{ 
	  ndeleted++;
	  continue;
	}

      (*itsex)->Iter()  = (*itdao)->Iter();
      (*itsex)->Chi()   = (*itdao)->Chi();
      (*itsex)->Sharp() = (*itdao)->Sharp();

      if (ReplacePos)
	{
	  (*itsex)->x = (*itdao)->x;
	  (*itsex)->y = (*itdao)->y;
	}

      if (ReplaceFlux)
	{
	  (*itsex)->flux = (*itdao)->flux;
	  (*itsex)->EFlux() = (*itdao)->EFlux();
	}

      if (ReplaceSky)(*itsex)->Fond() = (*itdao)->Fond();
    }   

  if (ndeleted>0) cout << " CombineSextractorDaophot : did not find "
		       << ndeleted << " in daophot catalog " << endl; 

  sexstars.write(DbIm.ImageCatalogName(SExtractor));
}


#endif // DAOPHOTIO__H
