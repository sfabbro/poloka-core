#include "senearstar.h"
#include "fitsimage.h"
#include "gtransfo.h"

//Implementation of class SENearStar


SENearStar::SENearStar(const SEStar &AStar, const int Dx, const int Dy, 
		       const double PhotFactor, const double GalFlux) 
  : SEStar(AStar)
{
  xShift =  Dx;
  yShift = Dy;
  photFactor = PhotFactor;
  galFlux = GalFlux;
} 

SENearStar::SENearStar() : SEStar()
{
  xShift =  0;
  yShift = 0;
  photFactor = 0;
  galFlux = 0;
}



SENearStar *SENearStar::ActualFake() const
{
  SENearStar *result = new SENearStar(*this);
  result->x += double(xShift);
  result->y += double(yShift);
  result->flux *= photFactor;
  return result;
}


string SENearStar::WriteHeader_(ostream &pr, const char *i) const
{
  if (i == NULL) i = "";
  string SEStarFormat = SEStar::WriteHeader_(pr,i);
  pr << "# dx"<< i << " : offset to real position" << endl;
  pr << "# dy"<< i << " : offset to real position" << endl;
  pr << "# photfac"<< i << " : phot factor from model to fake. " << endl;
  pr << "# galflux"<< i << " : galaxy flux on reference image " << endl;
  return SEStarFormat +" SENearStar 1 ";
}




void SENearStar::writen(ostream & pr) const
{
  SEStar::writen(pr);
  pr << ' ' << xShift << ' ' << yShift << ' ' <<  photFactor << ' ' << galFlux << ' ';
  //  return 1;
}

SENearStar* SENearStar::read(istream& r, const char* Format)
{
  SEStar *se = SEStar::read(r, Format);
  if (!se) return NULL;
  SENearStar *s = new SENearStar(*se);
  delete se;
  r >>  s->xShift >> s->yShift >> s->photFactor >> s->galFlux;
  return s;
}



#include "vutils.h"

static float sky_level(Image &image, Point &where, int stampSize)
{
  int xstart = int(where.x - stampSize);
  int ystart = int(where.y - stampSize);
  int xend = xstart + stampSize;
  int yend = ystart + stampSize;
  Pixel *values = new Pixel [4*stampSize*stampSize];
  int count = 0;
  int stampSize2 = stampSize * stampSize;
  for ( int j = ystart ; j < yend; ++j)
    {
      if (j < 0) continue;
      if (j >= image.Ny()) break;
      for ( int i = xstart ; i < xend; ++i)
      {
	if (i < 0) continue;
	if (i >= image.Ny()) break;
	  if (i*i + j*j < stampSize2) continue;
	  values[count] = image(i,j);
	  count ++;
	}
    }
  float median = FArrayMedian(values,count);
  delete [] values;
  return median;
}

   
void SENearStar::AddToFitsImage(Image &image, 
				Image & dest, const Gtransfo *Transfo) const
{
  // ou se trouve l'etoile modele.
  Point model = Transfo->apply(*this);	
  

  int stampSize = int(3.* Fwhm()+.5);

  // sky level around the model
  float skylevel = sky_level(image, model, 2*stampSize);

  //Determination des coordonnees de la vignette
  int xstart = int(model.x - stampSize + 0.5);
  int ystart = int(model.y - stampSize + 0.5);
  int xend = min(xstart + 2*stampSize,image.Nx());
  int yend = min(ystart + 2*stampSize,image.Ny());
  xstart = max(xstart,0);
  ystart = max(ystart,0);
  for (int j = ystart; j < yend; ++j)
    {
      int targetj = j+ yShift;
      if (targetj < 0) continue;
      if (targetj >= dest.Ny()) break;
      for (int i = xstart; i < xend; ++i)
	{
          int targeti = i + xShift;
	  if (targeti < 0 ) continue;
          if (targeti >= dest.Nx()) break;
	  dest(targeti , targetj) += (image(i,j)- skylevel)*photFactor;
	}
    }
  
} 



void SENearStarList::AddToFitsImage(FitsImage &image, const Gtransfo *Transfo) const
{
  for (SENearStarCIterator i = begin(); i != end(); ++i)
    {
      const SENearStar *s = *i;
      s->AddToFitsImage(image, Transfo);
    }
}


void SENearStarList::ActualFakes(SENearStarList &Result) const
{
  Result.clear();
  for (SENearStarCIterator it = begin();  it != end(); ++it)
      Result.push_back((*it)->ActualFake());
  }
// Converter :
BaseStarList* SENear2Base(SENearStarList * This)
{ return (BaseStarList*) This;}

const BaseStarList* SENear2Base(const SENearStarList * This)
{ return (BaseStarList*) This;} 

#include "starlist.cc" /* since starlist is a template class */
template class StarList<SENearStar>;  /* to force instanciation */
