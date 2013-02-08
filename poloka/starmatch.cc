#include <iostream>
#include <fstream>
#include <iomanip>

#include <poloka/gtransfo.h>
#include <poloka/starmatch.h>
#include <poloka/basestar.h>
#include <poloka/vutils.h>

/* TO DO:
   think about imposing a maximum number of matches that may 
   be discarded in Cleanup. 
*/

static double sq(const double &x) { return x*x;}

double StarMatch::Chi2(const Gtransfo &T) const
{
  FatPoint tr;
  T.TransformPosAndErrors(point1, tr);
  double vxx = tr.vx + point2.vx;
  double vyy = tr.vy + point2.vy;
  double vxy = tr.vxy + point2.vxy;
  double det = vxx*vyy-vxy*vxy;
  return (vyy*sq(tr.x-point2.x) + vxx*sq(tr.y-point2.y)
	  -2*vxy*(tr.x-point2.x)*(tr.y-point2.y))/det;
}


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

static double* chi2_array(const StarMatchList &L,const Gtransfo &T)
{
  unsigned s = L.size();
  double *res = new double[s];
  unsigned count = 0;
  for (StarMatchCIterator it= L.begin(); it!= L.end(); ++it)
    res[count++] = it->Chi2(T);
  return res;
}


static unsigned chi2_cleanup(StarMatchList &L,  const double Chi2Cut, 
			     const Gtransfo &T)
{
  unsigned erased = L.RemoveAmbiguities(T);
  for (StarMatchIterator smi = L.begin(); smi != L.end(); )
    {
      if (smi->chi2  > Chi2Cut)
	{
	  smi = L.erase(smi);
	  erased++;
	}
      else ++smi;
    }
  return erased;
}
  


static double *get_dist2_array(const StarMatchList &L, const Gtransfo &T)
{
  unsigned npair = L.size();
  if (npair == 0) return NULL; 
  double *dist = new double [npair];

  unsigned i=0;
  for (  StarMatchCIterator smi = L.begin(); smi != L.end(); ++smi, ++i)
    {
      const Point &p1 = smi->point1;
      const Point &p2 = smi->point2;
      dist[i] = p2.Dist2(T.apply(p1));
    }
  return dist;
}

int StarMatchList::Dof(const Gtransfo* T) const
{
  if (!T) T=transfo;
  int dof = 2*size() - T->Npar();
  return (dof>0) ? dof : 0;
}




  /*! removes pairs beyond NSigmas in distance (where the sigma scale is
     set by the fit) and iterates until stabilization of the number of pairs. 
     If the transfo is not assigned, it will be set to a GtransfoLinear. User
     can set an other type/degree using SetTransfo() before call. */
void StarMatchList::RefineTransfo(const double &NSigmas)
{
  double cut;
  unsigned nremoved;
  if (!transfo) transfo = new GtransfoLin;
  do 
    {
      int nused = size();
      if (nused <= 2) { chi2 = -1; break;}
      chi2 = transfo->fit(*this);
      /* convention of the fitted routines : 
	 -  chi2 = 0 means zero degrees of freedom 
              (this was not enforced in Gtransfo{Lin,Quad,Cub} ...)
	 -  chi2 = -1 means ndof <0 (and hence no possible fit)
	 --> in either case, refinement is over
         The fact that chi2 = 0 was not enforced when necessary means
	 that in this (rare) case, we were discarding matches at random....
         With GtransfoPoly::fit, this is no longer the case.
      */
      if (chi2 <= 0)  return;
      unsigned npair = int(size());
      if (npair == 0) break; // should never happen

      // compute some chi2 statistics
      double *chi2_array = new double[npair];
      unsigned count = 0;
      for (StarMatchIterator it= begin(); it!= end(); ++it)
	chi2_array[count++] = it->chi2 = it->Chi2(*transfo);
      double median = DArrayMedian(chi2_array,npair);
      delete [] chi2_array;

      // discard outliers : the cut is understood as a "distance" cut
      cut = sq(NSigmas)*median;
      nremoved = chi2_cleanup(*this, cut, *transfo);
    } 
  while (nremoved);
  dist2 = ComputeDist2(*this, *transfo);
}


/* not very robust : assumes that we went through Refine just before... */
double StarMatchList::Residual() const
{
  double deno = (2.*size() - transfo->Npar());
  return (deno>0.) ? sqrt(dist2/deno) : -1; // is -1 a good idea?
}


void StarMatchList::SetDistance(const Gtransfo &Transfo)
{
  for (StarMatchIterator smi = begin(); smi != end(); smi++) (*smi).SetDistance(Transfo); // c'est compact
}


unsigned StarMatchList::RemoveAmbiguities(const Gtransfo &Transfo, 
				     const int Which)
{
  if (!Which) return 0;
  SetDistance(Transfo);
  int initial_count = size();
  if (Which & 1)
    {
      sort(CompareS1);
      unique(SameS1);
    }
  if (Which & 2)
    {
      sort(CompareS2); unique(SameS2);
    }
  return (initial_count-size());
}



unsigned StarMatchList::Cleanup(double DistanceCut, const Gtransfo &ResultTransfo)
{
StarMatchIterator smi;

int erased = RemoveAmbiguities(ResultTransfo);
for (smi = begin(); smi != end(); )
  {
  double distance = (*smi).point2.Distance(ResultTransfo.apply((*smi).point1));
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
    case 2 : SetTransfo(new GtransfoPoly(2));     break;
    case 3 : SetTransfo(new GtransfoPoly(3));      break;
    case 4 : SetTransfo(new GtransfoPoly(4));      break;
    case 5 : SetTransfo(new GtransfoPoly(5));      break;
    case 6 : SetTransfo(new GtransfoPoly(6));      break;
    case 7 : SetTransfo(new GtransfoPoly(7));      break;
    case 8 : SetTransfo(new GtransfoPoly(8));      break;
    case 9 : SetTransfo(new GtransfoPoly(9));      break;
    case 10 : SetTransfo(new GtransfoPoly(10));      break;
    default : cerr << "Wrong transfo order : " << Order << endl; order = -1; return;
    }
  order = Order;
}



/* This routine should operate on a copy : RefineTransfo 
   might shorten the list */
Gtransfo* StarMatchList::InverseTransfo() /* it is not const although it tries not to change anything  */
{
  if (!transfo) return NULL;

  Gtransfo *old_transfo = transfo->Clone();  
  double old_chi2 = chi2;

  Swap();
  SetTransfoOrder(order);
  RefineTransfo(3.);// keep same order
  Gtransfo *inverted_transfo = transfo->Clone();
  SetTransfo(old_transfo);
  delete old_transfo;
  Swap();
  chi2 = old_chi2;

  return inverted_transfo;
}

 
void StarMatchList::CutTail(const int NKeep)
{
iterator si;
int count=0;
for (si = begin(); si != end() && count < NKeep; ++count, ++si);
erase(si, end());
}

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
      FatPoint p1 = *starm.s1;
      if (Transfo) 
	{
	  Transfo->TransformPosAndErrors(p1,p1);
	  pr << p1.x << ' ' << p1.y << ' ';
	  double sx = sqrt(p1.vx);
	  double sy = sqrt(p1.vy);
	  pr << sx << ' ' << sy << ' ' << p1.vxy/(sx*sy) << ' ';
	}
      (starm.s2)->writen(pr);

      // compute offsets here  because they can be rounded off by paw.
      double dx = p1.x - starm.s2->x;
      double dy = p1.y - starm.s2->y;
      pr << dx << ' '  << dy << ' ' << sqrt(dx*dx+dy*dy) << ' ';
      // chi2 assoc
      if (Transfo)
	pr << it->Chi2(*Transfo) << ' ';
      else 
	pr << it->Chi2(GtransfoIdentity()) << ' ';
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
  if (tf)
    {
      pr << "# x1tf: transformed x1 coordinate "  << endl ; 
      pr << "# y1tf: transformed y1 coordinate "  << endl ; 
      pr << "# sx1tf: transformed sx1 "  << endl ; 
      pr << "# sy1tf: transformed sy1 "  << endl ; 
      pr << "# rxy1tf: transformed rhoxy1 "  << endl ; 
    }
  (starm.s2)->WriteHeader_(pr, "2");
  pr << "# dx : diff in x" << std::endl;
  pr << "# dy : diff in y" << std::endl;
  pr << "# dass : association distance"  << endl ; 
  pr << "# chi2ass : assoc chi2"  << endl ; 
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
  for (StarMatchIterator i= begin(); i != end(); ++i)
    {
      i->chi2 = i->Chi2(*transfo);
      chi2 += i->chi2;
    }
}

void StarMatchList::DumpTransfo(ostream &stream) const
{ 
  stream << " ================================================================" << endl
	 << " Transformation between lists of order " << TransfoOrder() << endl
	 << *transfo //<< endl
	 << " Chi2 = " << Chi2() << "  Residual = " << Residual() << endl
	 << "  Number in the list = " << size() <<endl
	 << " ================================================================" << endl;
}


double FitResidual(const double Dist2, const StarMatchList &S, const Gtransfo &T)
{
  return sqrt(Dist2/(2.*S.size()-T.Npar()));
}


double FitResidual(const StarMatchList &S, const Gtransfo &T)
{
  return FitResidual(ComputeDist2(S,T),S,T);
}

double ComputeDist2(const StarMatchList &S, const Gtransfo &T)
{
  double dist2 = 0;
  for (StarMatchCIterator i = S.begin(); i != S.end(); ++i)
    dist2 += T.apply(i->point1).Dist2(i->point2);
  return dist2;
}

double ComputeChi2(const StarMatchList &L, const Gtransfo &T)
{
  unsigned s= L.size();
  double *chi2s = chi2_array(L,T);
  double chi2 = 0;
  for (unsigned k=0; k<s; ++k) chi2 += chi2s[k];
  delete [] chi2s;
  return chi2;
}
