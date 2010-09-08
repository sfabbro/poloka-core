#ifndef FAKELIST__CC
#define FAKELIST__CC
#include <iostream>
#include <iomanip>
#include <fstream>
#include <string>
#include <cstdio>
#include <cmath>

#include "fakelist.h"
#include "frame.h"

 

template<class Fake> FakeList<Fake>::FakeList(const string &FileName) /* to be changed if switch to Fake rather than pointers to Fakes */ 
{
  this->read(FileName);
}


template<class Fake> int FakeList<Fake>::read(istream & r)
{
  char car ;
  char buff[400];
  char *format = 0;
  ClearList();
  while( r >> car ) // to test eof
    {
      r.unget() ;
      if ( (car == '#') ) // we jump over the line  (not always ...)
        {
	  r.getline(buff,400);
	  /* hack something reading " format <FakeType> <integer>" to drive the decoding (in Fake::read) */
	  char *p = strstr(buff,"format");
	  if (p) /* this test is enough because the format is the last line of the header ... */
          {
	    format = p + strlen("format");
          }
        }
      else
	{
	  Fake* s = Fake::read(r, format); 
	  if (!s) 
	    {
	      return 0;
	    }
	  push_back(s);
	}
    }
  return 1 ;
}

template<class Fake> int 
FakeList<Fake>::read(const string &FileName)
{
  return ascii_read(FileName);
}

template<class Fake> int 
FakeList<Fake>::ascii_read(const string &FileName)
{
  ifstream rd(FileName.c_str());
  if (!rd)
    {
      cout << "FakeList cannot open :" << FileName << endl;
      return 0;
    }
  return read(rd);
}

template<class Fake> int FakeList<Fake>::write(const string &FileName) /* const */
{
  ofstream pr(FileName.c_str());
  if (!pr)
    {
      cerr << " FakeList::write : could not open " << FileName << endl;
      return 0;
    }
  write(pr);
  pr.close();
 return 1;
}


// check pr a faire avant, close a faire tout seul
template<class Fake> int FakeList<Fake>::write(ostream & pr) const
{
  ios::fmtflags  old_flags =  pr.flags(); 
  pr  << resetiosflags(ios::scientific) ;
  //pr  << setiosflags(0) ;
  pr  << setiosflags(ios::fixed) ;
  int oldprec = pr.precision();
  pr<< setprecision(8);
  const Fake *theFirst = *(this->begin());
  // if (!theFirst) return 0;
  if (this->size() == 0) // it seems that we can have this->size(0 == 0 with this->front() != NULL !
    {
      Fake dummy;
      dummy.WriteHeader(pr);
    }
  else theFirst->WriteHeader(pr);
  for (FakeCIterator it= this->begin(); it!=this->end() ; it++ )
    {    
      (*it)->write(pr);
    }
  pr.flags(old_flags);
  pr << setprecision(oldprec);
 return 1;
}

template<class Fake> int FakeList<Fake>::write_wnoheader(ostream & pr) const
{
  // check pr a faire avant, close a faire tout seul  
  ios::fmtflags  old_flags =  pr.flags(); 
  pr  << resetiosflags(ios::scientific) ;
  pr  << setiosflags(ios::fixed) ;
  int oldprec = pr.precision();
  pr<< setprecision(8);
  for (FakeCIterator it= this->begin(); it!= this->end() ; it++ )
    {    
      (*it)->Fake::write(pr);
    }
  pr.flags(old_flags);
  pr << setprecision(oldprec);
 return 1;
}






template <class Fake>  void FakeList<Fake>::ExtractHead(FakeList<Fake> &Out, int NHead) const
{
int count = 0;
for (FakeCIterator s= this->begin(); s!= this->end() && (count < NHead); ++s, count++)
  {
  Fake *copy = new Fake(*(*s));  /* to be changed if switch to Fake rather than pointers to Fakes */
  Out.push_back(copy);
  }
}



template<class Fake> Fake* FakeList<Fake>::FindClosest(double X, double Y) const
{
double min_dist2 = 1e30;
const Fake *minstar = NULL;
double dist2;
for (FakeCIterator si = this->begin(); si!= this->end(); ++si) 
   { 
   const Fake *s = *si;
   dist2 = (X - s->x)*(X - s->x) + (Y - s->y)*(Y - s->y);
   if (dist2 < min_dist2) { min_dist2 = dist2; minstar = s; }
   }
 return (Fake *) minstar; // violates constness
}

template<class Fake> Fake* FakeList<Fake>::FindClosest(double X, double Y, FakeList<Fake> const & stlparallel, Fake* & star2  ) const
{
  if (this->size() != stlparallel.size() )
    cerr << "Size list differents in FindClosest" << endl ;
double sqmin_dist = 1e30;
const Fake *minstar = NULL;
double sqdist;
FakeCIterator si2 = stlparallel.begin(); 
const Fake *minstar2 = NULL ;
for (FakeCIterator si = this->begin(); si!= this->end() && si2 != stlparallel.end(); ++si) 
   { 
   const Fake *s = *si;
   const Fake *s2 = *si2;
   sqdist = (X - s->x)*(X - s->x) + (Y - s->y)*(Y - s->y);
   if (sqdist < sqmin_dist) 
     { 
       sqmin_dist = sqdist; 
       minstar = s; 
       minstar2 = s2 ; 
     }
   ++si2 ;
   }
 star2 = (Fake *) minstar2; 
 return (Fake *) minstar; // violates constness
}


template<class Fake> Fake* FakeList<Fake>::FindAndRemoveClosest(double X, double Y) 
{
double min_dist2 = 1e30;
const Fake *minstar = NULL;
double dist2;
 FakeIterator si_res;
for (FakeIterator si = this->begin(); si!= this->end(); ++si) 
   { 
   const Fake *s = *si;
   dist2 = (X - s->x)*(X - s->x) + (Y - s->y)*(Y - s->y);
   if (dist2 < min_dist2) { min_dist2 = dist2; minstar = s;
   si_res = si ;}
   }
 if ( minstar == NULL)
 return (Fake *) minstar; 
 else
   {
     Fake * star_res = new Fake(*minstar);
     this->erase(si_res);
     return (Fake *) star_res; // violates constness
   }
}

template<class Fake> bool FakeList<Fake>::HasCloseNeighbor(double X, double Y, double maxdist, double mindist) const
{
  double dist2;
  double mindist2 = mindist*mindist;
  double maxdist2 = maxdist*maxdist;
  for (FakeCIterator si = this->begin(); si!= this->end(); ++si) 
    { 
      const Fake *s = (*si); 
      dist2 = (X - s->x)*(X - s->x) + (Y - s->y)*(Y - s->y);
      bool okflux = true ;
      if (okflux && (dist2 > mindist2) && (dist2 < maxdist2)) {
	cerr << "HasCloseNeighbor at dist : " << sqrt(dist2) << " Neighbor: " ;
	s->dump();
	return true;}
    }
  return false;
}  

template<class Fake> Fake* FakeList<Fake>::ClosestNeighbor(double X, double Y, double mindist) const
{
  const Fake *minstar = NULL;
  double dist2;
  double mindist2 = mindist*mindist;
  double min_dist2 = 1e30;
  for (FakeCIterator si = this->begin(); si!= this->end(); ++si) 
    { 
      const Fake *s = (*si); 
      dist2 = (X - s->x)*(X - s->x) + (Y - s->y)*(Y - s->y);
      if ((dist2 > mindist2) && (dist2 < min_dist2)) { min_dist2 = dist2; minstar = s;}
    }
  return (Fake *) minstar;  // violates constness
}


template<class Fake> int FakeList<Fake>::NumberOfNeighbors(const double &X, const double &Y, const double &distmax) const
{
  int nstars = 0;
  double dist2 = distmax*distmax;

  for (FakeCIterator it = this->begin(); it!= this->end(); ++it) 
    { 
      const Fake *s  = *it; 
      if ((X - s->x)*(X - s->x) + (Y - s->y)*(Y - s->y) < dist2) 
	{
	  ++nstars;
	}
    }
  return nstars;
}

template<class Fake> int FakeList<Fake>::AllNeighbors(FakeList &NeighborList, const double &X, 
						      const double &Y, const double &distmax) const
{
  int nstars = 0;
  double dist2 = distmax*distmax;
  NeighborList.ClearList();
  for (FakeCIterator it = this->begin(); it!= this->end(); ++it) 
    { 
      const Fake *s  = *it; 
      if ((X - s->x)*(X - s->x) + (Y - s->y)*(Y - s->y) < dist2) 
	{
	  NeighborList.push_back(new Fake(*s));
	  ++nstars;
	}
    }
  return nstars;
}


template<class Fake>void FakeList<Fake>::CutTail(const int NKeep)
{
int count = 0;
FakeIterator si;
for (si = this->begin(); si != this->end() && count < NKeep; ++count, ++si);
while ( si != this->end() ) {si = erase(si);}
}


template<class Fake>void FakeList<Fake>::ExtractInFrame(FakeList<Fake> &Out, const Frame &aFrame) const
{
  for (FakeCIterator s= this->begin(); s!= this->end(); ++s)
    {
      const Fake *st  = *s;
      if (aFrame.InFrame(st->x,st->y))
	{
	  Fake *copy = new Fake(*st);
	  Out.push_back(copy);
	}
    }
}

template<class Fake>void FakeList<Fake>::CutEdges(const Frame &aFrame, float mindist) 
{
  for (FakeIterator si= this->begin(); si!= this->end();)
    {
      if (aFrame.MinDistToEdges(**si) < mindist)
	{
	  si = erase(si); 
	}
      else
	si++;
    }
}

template<class Fake>void FakeList<Fake>::CopyTo(FakeList<Fake> &Copy) const
{
  Copy.ClearList();
  FakeCIterator si;
  for (si = this->begin(); si != this->end(); ++si) Copy.push_back(new Fake(*(*si)));
}

#ifdef STORAGE
template<class Fake> Fake*  FakeList<Fake>::CloneFirst() const
{
Fake *star = *(this->begin());
 return ( (star)? new Fake(*star) : NULL);
}
#endif

#ifdef USE_ROOT
#include "rootstuff.h"
#include <TClass.h>
/* old fashioned home made Streamer, because rootcint does not handle
the inheritance scheme properly ( FakeListWithRoot<> : FakeList<> : list<> )
rootcint does not realize that FakeListWithRoot is a list...
*/

template<class Fake> void FakeListWithRoot<Fake>::Streamer(TBuffer &R__b)
{
   // Stream an object of class 

   if (R__b.IsReading()) {
      FakeListWithRoot<Fake>::Class()->ReadBuffer(R__b, this);
      {
         ClearList();
         int R__i, R__n;
         R__b >> R__n;
         for (R__i = 0; R__i < R__n; R__i++) {
            Fake* R__t;
            R__b >> R__t;
            push_back(R__t);
         }
      }
   } else {
      FakeListWithRoot<Fake>::Class()->WriteBuffer(R__b, this);
      {
         R__b << int(size());
         FakeList<Fake>::FakeCIterator R__k;
         for (R__k = begin(); R__k != end(); ++R__k) {
	   const Fake* s = *R__k;
            R__b << s;
         }
      }
   }
}

template <class Fake> FakeListWithRoot<Fake>::FakeListWithRoot() : FakeList<Fake>() 
{};

template <class Fake> FakeListWithRoot<Fake>::FakeListWithRoot(const string &FileName)
{
  this->read(FileName);
}

template <class Fake> int FakeListWithRoot<Fake>::read(const string &FileName)
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
template <class Fake> int FakeListWithRoot<Fake>::write(const string& FileName)
{
  return root_write(FileName);
}



#include <TFile.h>
#include <TKey.h> 

template <class Fake> int FakeListWithRoot<Fake>::root_read(const string &FileName)
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
template <class Fake> int FakeListWithRoot<Fake>::root_write(const string &FileName)
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
#endif /* FAKELIST__CC */
