#ifndef STARLIST__CC
#define STARLIST__CC
#include <iostream>
#include <iomanip>
#include <fstream>
#include <string>
#include <cstdio>


#include "starlist.h"
#include "frame.h"



template<class Star> StarList<Star>::StarList(const string &FileName) /* to be changed if switch to Star rather than pointers to Stars */ 
{
  this->read(FileName);
}


template<class Star> int StarList<Star>::read(istream & r)
{
  char c ;
  char buff[400];
  char *format = 0;
  ClearList();
  while( r >> c ) // to test eof
    {
      r.unget() ;
      if ( (c == '#') ) // we jump over the line  (not always ...)
        {
	  r.getline(buff,400);
	  /* hack something reading " format <StarType> <integer>" to drive the decoding (in Star::read) */
	  char *p = strstr(buff,"format");
	  if (p) /* this test is enough because the format is the last line of the header ... */
          {
	    format = p + strlen("format");
          }
        }
      else
	{
	  Star* s = Star::read(r, format); 
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
  ifstream rd(FileName.c_str());
  if (!rd)
    {
      cout << "StarList cannot open :" << FileName << endl;
      return 0;
    }
  return read(rd);
}

template<class Star> int StarList<Star>::write(const string &FileName) /* const */
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
  // check pr a faire avant, close a faire tout seul
  const Star *theFirst = *(begin());
  // if (!theFirst) return 0;
  if (this->size() == 0) // it seems that we can have this->size(0 == 0 with this->front() != NULL !
    {
      Star dummy;
      dummy.WriteHeader(pr);
    }
  else theFirst->WriteHeader(pr);
  for (StarCIterator it= begin(); it!=end() ; it++ )
    {    
      (*it)->write(pr);
    }
  pr.flags(old_flags);
  pr << setprecision(oldprec);
 return 1;
}

template<class Star> int StarList<Star>::write_wnoheader(ostream & pr) const
{
  // check pr a faire avant, close a faire tout seul
  for (StarCIterator it= begin(); it!=end() ; it++ )
    {    
      (*it)->Star::write(pr);
    }
 return 1;
}






template <class Star>  void StarList<Star>::ExtractHead(StarList<Star> &Out, int NHead) const
{
int count = 0;
for (StarCIterator s= begin(); s!=end() && (count < NHead); ++s, count++)
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
for (StarCIterator si = begin(); si!= end(); ++si) 
   { 
   const Star *s = *si;
   dist2 = (X - s->x)*(X - s->x) + (Y - s->y)*(Y - s->y);
   if (dist2 < min_dist2) { min_dist2 = dist2; minstar = s;}
   }
 return (Star *) minstar; // violates constness
}

template<class Star>bool StarList<Star>::HasCloseNeighbor(double X, double Y, double maxdist, double mindist) const
{
  double dist2;
  double mindist2 = mindist*mindist;
  double maxdist2 = maxdist*maxdist;
  for (StarCIterator si = begin(); si!= end(); ++si) 
    { 
      const Star *s = (*si); 
      dist2 = (X - s->x)*(X - s->x) + (Y - s->y)*(Y - s->y);
      
      if ((dist2 > mindist2) && (dist2 < maxdist2)) {return true;}
    }
  return false;
}  

template<class Star> Star* StarList<Star>::ClosestNeighbor(double X, double Y, double mindist) const
{
  const Star *minstar = NULL;
  double dist2;
  double mindist2 = mindist*mindist;
  double min_dist2 = 1e30;
  for (StarCIterator si = begin(); si!= end(); ++si) 
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

  for (StarCIterator it = begin(); it!= end(); ++it) 
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
  for (StarCIterator it = begin(); it!= end(); ++it) 
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


template<class Star>void StarList<Star>::CutTail(const int NKeep)
{
int count = 0;
StarIterator si;
for (si = begin(); si != end() && count < NKeep; ++count, ++si);
while ( si !=end() ) {si = erase(si);}
}


template<class Star>void StarList<Star>::ExtractInFrame(StarList<Star> &Out, const Frame &aFrame) const
{
  for (StarCIterator s= begin(); s!=end(); ++s)
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
  for (StarIterator si= begin(); si!=end();)
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
  for (si = begin(); si != end(); ++si) Copy.push_back(new Star(*(*si)));
}

#ifdef STORAGE
template<class Star> Star*  StarList<Star>::CloneFirst() const
{
Star *star = *(this->begin());
 return ( (star)? new Star(*star) : NULL);
}
#endif

#ifdef USE_ROOT
#include "rootstuff.h"
#include <TClass.h>
/* old fashioned home made Streamer, because rootcint does not handle
the inheritance scheme properly ( StarListWithRoot<> : StarList<> : list<> )
rootcint does not realize that StarListWithRoot is a list...
*/

template<class Star> void StarListWithRoot<Star>::Streamer(TBuffer &R__b)
{
   // Stream an object of class 

   if (R__b.IsReading()) {
      StarListWithRoot<Star>::Class()->ReadBuffer(R__b, this);
      {
         ClearList();
         int R__i, R__n;
         R__b >> R__n;
         for (R__i = 0; R__i < R__n; R__i++) {
            Star* R__t;
            R__b >> R__t;
            push_back(R__t);
         }
      }
   } else {
      StarListWithRoot<Star>::Class()->WriteBuffer(R__b, this);
      {
         R__b << int(size());
         StarList<Star>::StarCIterator R__k;
         for (R__k = begin(); R__k != end(); ++R__k) {
	   const Star* s = *R__k;
            R__b << s;
         }
      }
   }
}

template <class Star> StarListWithRoot<Star>::StarListWithRoot() : StarList<Star>() 
{};

template <class Star> StarListWithRoot<Star>::StarListWithRoot(const string &FileName)
{
  this->read(FileName);
}

template <class Star> int StarListWithRoot<Star>::read(const string &FileName)
{
  FILE *file = fopen(FileName.c_str(), "r");
  if (!file) 
    {
      cerr << " Cannot open catalog " << FileName << endl;
      return 0;
    }
  char start[4];
  if (fread(start,1,1,file) != 1)
    {
      cerr << " empty calalog file " << FileName << endl;
      fclose(file);
      return 0;
    }
  fclose(file);
  if (start[0] == '#') ascii_read(FileName);
  else root_read(FileName);
  return this->size();
}

/* we may just rename root_write */
template <class Star> int StarListWithRoot<Star>::write(const string& FileName)
{
  return root_write(FileName);
}



#include <TFile.h>
#include <TKey.h> 

template <class Star> int StarListWithRoot<Star>::root_read(const string &FileName)
{
  TFile tfile(FileName.c_str());
  TIter nextkey(tfile.GetListOfKeys());
  TKey *key = (TKey*)nextkey();
  if (key) 
    {
      key->Read(this);
    }
  else 
    {
      cerr << " could not read : " << FileName << endl;
    }
  tfile.Close(); 
  return this->size();      
}

//cannot be const since TObject::Write is not
template <class Star> int StarListWithRoot<Star>::root_write(const string &FileName)
{
  TFile tfile(FileName.c_str(),"RECREATE");
  if (!tfile.IsOpen())
    {
      cerr << " could not open " << FileName << endl;
      return 0;
    }
  int count = Write();
  tfile.Close();
  return count;
}
#endif /* USE_ROOT */
#endif /* STARLIST__CC */
