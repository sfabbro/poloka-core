#include <iostream>
#include <fstream>
#include <iomanip>

#include "gtransfo.h"
#include "starmatch.h"
#include "basestar.h"
#include "algorithm" // for copy
#include "vutils.h" /* for DArrayMedian */

ostream& operator << (ostream &stream, const StarMatch &Match)
{ 
  stream << Match.point1.x << ' ' << Match.point1.y << ' ' 
	 << Match.point2.x << ' ' << Match.point2.y << ' ' 
	 << Match.distance << endl; return stream; 
}

ostream& operator << (ostream &stream, const StarMatchList &List)
{ 
  stream << " number of elements " << List.size() << endl; 
  copy(List.begin(), List.end(), ostream_iterator<StarMatch>(stream)); 
  return stream;
} 

double *StarMatchList::Dist2() const
{
  StarMatchCIterator smi;
  int npair = int(this->size());
  if (npair == 0) return NULL; 
  double *dist = new double [npair];

  int i=0;
  for (smi = begin(); smi != end(); ++smi, ++i)
    {
      const StarMatch &match = *smi;
      const Point &p1 = match.point1;
      const Point &p2 = match.point2;
      dist[i] = p2.Dist2(transfo->apply(p1));
    }
  return dist;
}

  /*! removes pairs beyond NSigmas (where the sigma scale is
     set by the fit) and iterates until stabilization of the number of pairs. 
     If the transfo is not assigned, it will be set to a GtransfoLinear. User
     can set an other type/degree using SetTransfo() before call. */
void StarMatchList::RefineTransfo(const double &NSigmas)
{
  double cut;
  int nremoved;
  if (!transfo) transfo = new GtransfoLin;
  do 
    {
      nused = size();
      if (nused <= 2) { chi2 = -1; break;}
      chi2 = transfo->fit(*this);
      // cerr << *transfo << endl; see cerr instead of cout : result will change if sextractor pass before!!
      if (chi2<0)
	{
	  // cerr << " in  StarMatchList::Refine , chi2 = " << chi2 << endl;
	  return;
	}
      int npair = int(size());
      if (npair == 0) break;
      double *dist = Dist2();
      double median = DArrayMedian(dist,npair);
      delete [] dist;
      cut = NSigmas*sqrt(median);
      //cout << " nused " << nused << " chi2 " << chi2 << endl;
      nremoved = Cleanup(cut, *transfo);
      //cout << *transfo << endl;
    } 
  while (nremoved);
}



void StarMatchList::SetDistance(const Gtransfo &Transfo)
{
  for (StarMatchIterator smi = begin(); smi != end(); smi++) (*smi).SetDistance(Transfo); // c'est compact
}


int StarMatchList::RemoveAmbiguities(const Gtransfo &Transfo)
{
  SetDistance(Transfo);
  int initial_count = size();
  sort(CompareS1);
  //cout << " apres compare \n " << *this << endl;
  unique(SameS1);
  //cout << " apres unique \n " << *this << endl;
  sort(CompareS2); unique(SameS2);
  return (initial_count-size());
}



int StarMatchList::Cleanup(double DistanceCut, const Gtransfo &ResultTransfo)
{
StarMatchIterator smi;

int erased = RemoveAmbiguities(ResultTransfo);
for (smi = begin(); smi != end(); )
  {
  double distance = (*smi).point2.Distance(ResultTransfo(((*smi).point1)));
  if (distance > DistanceCut)
    {
    smi = erase(smi);
    erased++;
    }
  else ++smi;
  }
return erased;
}

void StarMatchList::SetTransfoOrder(const int Order)
{
  switch (Order)
    {
    case 0 : SetTransfo(new GtransfoLinShift()); break;
    case 1 : SetTransfo(new GtransfoLin());      break;
    case 2 : SetTransfo(new GtransfoQuad());     break;
    case 3 : SetTransfo(new GtransfoCub());      break;
    default : cerr << "Wrong transfo order : " << Order << endl; order = -1; return;
    }
  order = Order;
}


Gtransfo* StarMatchList::InverseTransfo() /* it is not const although it tries not to change anything  */
{
  if (!transfo) return NULL;

  Gtransfo *old_transfo = transfo->Clone();  
  double old_chi2 = chi2;
  int old_nused = nused;

  Swap();
  SetTransfoOrder(order);
  RefineTransfo(3.);// keep same order
  Gtransfo *inverted_transfo = transfo->Clone();
  SetTransfo(old_transfo);
  delete old_transfo;
  Swap();
  chi2 = old_chi2;
  nused = old_nused;

  return inverted_transfo;
}

 
void StarMatchList::CutTail(const int NKeep)
{
iterator si;
int count=0;
for (si = begin(); si != end() && count < NKeep; ++count, ++si);
erase(si, end());
}

#ifdef STORAGE
void StarMatchList::write(ostream & pr) const
{

  if ( empty() )
    {
      cerr << " Can't write empty StarMatchList " << endl ;
      return ;
    }

  const StarMatch &starm = front(); 
  (starm.s1)->WriteHeader_(pr,"1");
  (starm.s2)->WriteHeader_(pr,"2");
  pr << "# dx : diff in x" << std::endl;
  pr << "# dy : diff in y" << std::endl;
  pr << "# dass : association distance"  << endl ; 
  pr << "# end " << endl ;

 
  GtransfoIdentity id ;
  for (StarMatchCIterator it= begin(); it!= end(); it++ )
    {
      const StarMatch &starm = *it ;

      (starm.s1)->writen(pr);
      pr << " " ;
      (starm.s2)->writen(pr);
      pr << " " ;
      double dx = starm.s1->x - starm.s2->x;
      double dy = starm.s1->y - starm.s2->y;
      pr << dx << ' '  << dy << ' ' << sqrt(dx*dx+dy*dy);
      pr << endl ;
    }
}
#endif


void StarMatchList::write_wnoheader(ostream & pr, const Gtransfo* Transfo) const
{

  if ( empty() )
    {
      cerr << " Can't write empty StarMatchList " << endl ;
      return ;
    }

  ios::fmtflags  old_flags =  pr.flags(); 
  pr  << resetiosflags(ios::scientific) ;
  pr  << setiosflags(ios::fixed) ;
  int oldprec = pr.precision();
  pr<< setprecision(10);
  for (StarMatchCIterator it= begin(); it!= end(); it++ )
    {
      StarMatch starm = *it ;

      (starm.s1)->writen(pr);
      pr << " " ;
      // transformed coordinates
      Point p1 = *starm.s1;
      if (Transfo) p1=Transfo->apply(p1);
      pr << p1.x << ' ' << p1.y << ' ';

      (starm.s2)->writen(pr);

      // compute offsets here  because they can be rounded off by paw.
      double dx = p1.x - starm.s2->x;
      double dy = p1.y - starm.s2->y;
      pr << dx << ' '  << dy << ' ' << sqrt(dx*dx+dy*dy);
      pr << endl ;
    }
  pr.flags(old_flags);
  pr << setprecision(oldprec);
}

void StarMatchList::write(ostream &pr, const Gtransfo *tf) const
{
  if ( empty() )
    {
      cerr << " Can't write empty StarMatchList " << endl ;
      return ;
    }

  const StarMatch &starm = front(); 
  (starm.s1)->WriteHeader_(pr, "1");
  pr << "# x1tf: coordonnees x1 transformee "  << endl ; 
  pr << "# y1tf: coordonnees y1 transformee "  << endl ; 

  (starm.s2)->WriteHeader_(pr, "2");
  pr << "# dx : diff in x" << std::endl;
  pr << "# dy : diff in y" << std::endl;
  pr << "# dass : association distance"  << endl ; 
  pr << "# end " << endl ;

  write_wnoheader(pr, tf);
}

void
StarMatchList::write(const string &filename, const Gtransfo *tf) const
{
  ofstream pr(filename.c_str()) ;
  write(pr, tf);
  pr.close();
}

#ifdef STORAGE
void
StarMatchList::write(const string &filename) const
{
  ofstream pr(filename.c_str()) ;
  //  cerr <<"Writing with fixed precision " << endl ;
  pr  << resetiosflags(ios::scientific) ;
  pr  << setiosflags(ios::fixed) ;
  int oldprec = pr.precision();
  pr<< setprecision(10);
  write(pr);
  pr.close();
}
#endif






void StarMatchList::Swap()
{
  for (StarMatchIterator it= begin(); it!= end(); ++it )
    {
      it->Swap() ; 
    }
}

int StarMatchList::RecoveredNumber(double mindist) const
{
  int n = 0 ;
  GtransfoIdentity identity;
  for (StarMatchCIterator it= begin(); it!= end(); ++it )
    {
      if ((*it).Distance(identity) < mindist)
	n++ ;
    }
  return(n);
}


void StarMatchList::ApplyTransfo(StarMatchList &Transformed, 
				 const Gtransfo *PriorTransfo,
				 const Gtransfo *PosteriorTransfo) const
{
  Transformed.clear();
  GtransfoIdentity id;
  const Gtransfo &T1 = (PriorTransfo)? *PriorTransfo : id;
  const Gtransfo &T2 = (PosteriorTransfo)? *PosteriorTransfo : id;

  for (StarMatchCIterator it= begin(); it!= end(); ++it )
    {
      Point p1 = T1.apply(it->point1);
      Point p2 = T2.apply(it->point2);
      Transformed.push_back(StarMatch(p1, p2, it->s1, it->s2));
    }
}

void StarMatchList::SetChi2()
{
  chi2 = 0;
  nused = int(size());
  double *dist2 = Dist2();
  for (int i=0; i<nused; ++i) chi2 += dist2[i];
  if (dist2) delete [] dist2;
}

void StarMatchList::DumpTransfo(ostream &stream) const
{ 
  stream << " ================================================================" << endl
	 << " Transformation between lists of order " << TransfoOrder() << endl
	 << *transfo //<< endl
	 << " Chi2 = " << Chi2() << "  Residual = " << Residual() << endl
	 << " Nused = " << Nused() <<  "  Number in the list = " << size() <<endl
	 << " ================================================================" << endl;
}


double FitResidual(const double Chi2, const StarMatchList &S, const Gtransfo &T)
{
  return sqrt(Chi2/(2.*S.size()-T.Npar()));
}



