#ifndef STARLIST__CC
#define STARLIST__CC
#include <iostream>
#include <iomanip>
#include <fstream>
#include <string>
#include <cstdio>


#include "starlist.h"
#include "frame.h"
#include "fastifstream.h"


template<class Star> StarList<Star>::StarList(const string &FileName) /* to be changed if switch to Star rather than pointers to Stars */ 
{
  this->read(FileName);
}


template<class Star> int StarList<Star>::read(fastifstream & r)
{
  char c ;
  char buff[4096];
  char format[512];
  ClearList();
  while( r >> c ) // to test eof
    {
      r.unget() ;
      if ( (c == '@') ) 
	{
	  r.getline(buff,4096); 
	  glob.ProcessLine(buff);
	  continue;
	}
      if ( (c == '#') ) // we jump over the line  (not always ...)
        {
	  r.getline(buff,400);
	  /* hack something reading " format <StarType> <integer>" to drive the decoding (in Star::read) */
	  char *p = strstr(buff,"format");
	  if (p) /* this test is enough because the format is the last line of the header ... */
          {
	    strncpy(format,p + strlen("format"),512);
          }
        }
      else
	{
	  Star* s = dynamic_cast<Star*>(Star::read(r, format));
	  // read the end of line, in case there are some items left
	  r.getline(buff,4096);
	  if (!s) 
	    {
	      return 0;
	    }
	  push_back(s);
	}
    }
  return 1 ;
}

template<class Star> int 
StarList<Star>::read(const string &FileName)
{
  return ascii_read(FileName);
}

template<class Star> int 
StarList<Star>::ascii_read(const string &FileName)
{
  fastifstream rd(FileName.c_str());
  if (!rd)
    {
      cout << "StarList cannot open :" << FileName << endl;
      return 0;
    }
  return read(rd);
}

template<class Star> int StarList<Star>::write(const string &FileName) const
{
  ofstream pr(FileName.c_str());
  if (!pr)
    {
      cerr << " StarList::write : could not open " << FileName << endl;
      return 0;
    }
  write(pr);
  pr.close();
 return 1;
}



template<class Star> int StarList<Star>::write(ostream & pr) const
{
  ios::fmtflags  old_flags =  pr.flags(); 
  pr  << resetiosflags(ios::scientific) ;
  // pr  << setiosflags(0) ;
  pr  << setiosflags(ios::fixed) ;
  int oldprec = pr.precision();
  pr<< setprecision(8);

  // write GlobalValues if any
  vector<string> globs = glob.OutputLines();
  for (unsigned k = 0; k < globs.size(); ++k)
    pr << '@' << globs[k] << endl;

  // check pr a faire avant, close a faire tout seul
  const Star *theFirst = this->front();
  // if (!theFirst) return 0;
  if (this->size() == 0) // it seems that we can have this->size(0 == 0 with this->front() != NULL !
    {
      Star dummy;
      dummy.WriteHeader(pr);
    }
  else theFirst->WriteHeader(pr);
  for (StarCIterator it= this->begin(); it!=this->end() ; it++ )
    {    
      (*it)->write(pr);
    }
  pr.flags(old_flags);
  pr << setprecision(oldprec);
 return 1;
}






template <class Star>  void StarList<Star>::ExtractHead(StarList<Star> &Out, int NHead) const
{
int count = 0;
for (StarCIterator s= this->begin(); s!= this->end() && (count < NHead); ++s, count++)
  {
  Star *copy = new Star(*(*s));  /* to be changed if switch to Star rather than pointers to Stars */
  Out.push_back(copy);
  }
}

template<class Star> Star* StarList<Star>::FindClosest(double X, double Y) const
{
double min_dist2 = 1e30;
const Star *minstar = NULL;
double dist2;
for (StarCIterator si = this->begin(); si!= this->end(); ++si) 
   { 
   const Star *s = *si;
   dist2 = (X - s->x)*(X - s->x) + (Y - s->y)*(Y - s->y);
   if (dist2 < min_dist2) { min_dist2 = dist2; minstar = s;}
   }
 return (Star *) minstar; // violates constness
}
template<class Star> Star* StarList<Star>::FindAndRemoveClosest(double X, double Y) 
{
double min_dist2 = 1e30;
const Star *minstar = NULL;
double dist2;
 StarIterator si_res;
for (StarIterator si = this->begin(); si!= this->end(); ++si) 
   { 
   const Star *s = *si;
   dist2 = (X - s->x)*(X - s->x) + (Y - s->y)*(Y - s->y);
   if (dist2 < min_dist2) { min_dist2 = dist2; minstar = s;
   si_res = si ;}
   }
 if ( minstar == NULL)
 return (Star *) minstar; 
 else
   {
     Star * star_res = new Star(*minstar);
     this->erase(si_res);
     return (Star *) star_res; // violates constness
   }
}

template<class Star>bool StarList<Star>::HasCloseNeighbor(double X, double Y, double maxdist, double mindist, double minflux) const
{
  double dist2;
  double mindist2 = mindist*mindist;
  double maxdist2 = maxdist*maxdist;
  for (StarCIterator si = this->begin(); si!= this->end(); ++si) 
    { 
      const Star *s = (*si); 
      dist2 = (X - s->x)*(X - s->x) + (Y - s->y)*(Y - s->y);
      bool okflux = true ;
      if (minflux > 0 )
	okflux = (s->flux > minflux);
      if (okflux && (dist2 > mindist2) && (dist2 < maxdist2)) {
	//cerr << "HasCloseNeighbor at dist : " << sqrt(dist2) << " Neighbor: " ;
	//s->dump();
	return true;}
    }
  return false;
}  

template<class Star> Star* StarList<Star>::ClosestNeighbor(double X, double Y, double mindist) const
{
  const Star *minstar = NULL;
  double dist2;
  double mindist2 = mindist*mindist;
  double min_dist2 = 1e30;
  for (StarCIterator si = this->begin(); si!= this->end(); ++si) 
    { 
      const Star *s = (*si); 
      dist2 = (X - s->x)*(X - s->x) + (Y - s->y)*(Y - s->y);
      if ((dist2 > mindist2) && (dist2 < min_dist2)) { min_dist2 = dist2; minstar = s;}
    }
  return (Star *) minstar;  // violates constness
}


template<class Star> int StarList<Star>::NumberOfNeighbors(const double &X, const double &Y, const double &distmax) const
{
  int nstars = 0;
  double dist2 = distmax*distmax;

  for (StarCIterator it = this->begin(); it!= this->end(); ++it) 
    { 
      const Star *s  = *it; 
      if ((X - s->x)*(X - s->x) + (Y - s->y)*(Y - s->y) < dist2) 
	{
	  ++nstars;
	}
    }
  return nstars;
}

template<class Star> int StarList<Star>::AllNeighbors(StarList &NeighborList, const double &X, 
						      const double &Y, const double &distmax) const
{
  int nstars = 0;
  double dist2 = distmax*distmax;
  NeighborList.ClearList();
  for (StarCIterator it = this->begin(); it!= this->end(); ++it) 
    { 
      const Star *s  = *it; 
      if ((X - s->x)*(X - s->x) + (Y - s->y)*(Y - s->y) < dist2) 
	{
	  NeighborList.push_back(new Star(*s));
	  ++nstars;
	}
    }
  return nstars;
}

template<class Star>
bool DecreasingFlux(const Star *S1, const Star *S2)
{
  return (S1->flux > S2->flux);
}

template<class Star>void StarList<Star>::FluxSort()
{
  this->sort(&DecreasingFlux<Star>);
}

template<class Star>void StarList<Star>::CutTail(const int NKeep)
{
int count = 0;
StarIterator si;
for (si = this->begin(); si != this->end() && count < NKeep; ++count, ++si);
while ( si != this->end() ) {si = erase(si);}
}


template<class Star>void StarList<Star>::ExtractInFrame(StarList<Star> &Out, const Frame &aFrame) const
{
  for (StarCIterator s= this->begin(); s!= this->end(); ++s)
    {
      const Star *st  = *s;
      if (aFrame.InFrame(*st))
	{
	  Star *copy = new Star(*st);
	  Out.push_back(copy);
	}
    }
}

template<class Star>void StarList<Star>::CutEdges(const Frame &aFrame, float mindist) 
{
  for (StarIterator si= this->begin(); si!= this->end();)
    {
      if (aFrame.MinDistToEdges(**si) < mindist)
	{
	  si = erase(si); 
	}
      else
	si++;
    }
}

template<class Star>void StarList<Star>::CopyTo(StarList<Star> &Copy) const
{
  Copy.ClearList();
  StarCIterator si;
  for (si = this->begin(); si != this->end(); ++si) Copy.push_back(new Star(*(*si)));
}

#ifdef STORAGE
template<class Star> Star*  StarList<Star>::CloneFirst() const
{
Star *star = *(this->begin());
 return ( (star)? new Star(*star) : NULL);
}
#endif


#endif /* STARLIST__CC */
