#include <iostream>
#include <fstream>

#include "gtransfo.h"
#include "nstarmatch.h"
#include "basestar.h"
#include "vutils.h" /* for DArrayMedian */
#include "fastfinder.h"

TStarList::TStarList(const BaseStarList& slist, const Gtransfo& transfo) {
  for(BaseStarCIterator it = slist.begin(); it != slist.end(); ++it) {
    push_back(new TStar(**it,transfo));
  }
}

string TStar::WriteHeader_(ostream & stream, const char*i) const
{ 
  return original->WriteHeader_(stream,i);
};


void TStar::writen(ostream &s) const 
{
  original->writen(s);
}


//==============================================================

void NStarMatch::AddMatch(const BaseStar &Match, const unsigned Index)
{
  // fill the array with NULLS up to Index-1
  // since Index is unsigned, do not compute Index-1 !!!!
  while (stars.size() + 1  <= Index )
    {
      stars.push_back(CountedRef<BaseStar>());
    }
  if (stars.size() > Index) // we may already have a match
    {
      DropMatch(Index);
      stars[Index] = CountedRef<BaseStar>(&Match);
    }
  else stars.push_back(CountedRef<BaseStar>(&Match));
  x = (nm*x+Match.x)/(nm+1);
  y = (nm*y+Match.y)/(nm+1);
  nm++;
}


void NStarMatch::DropMatch(const unsigned Index)
{
  if (Index >= stars.size()) return;
  BaseStarRef &toDrop = stars[Index];
  if (toDrop == NULL) return;
  // update position 
  if (nm > 1)
    {
      x = (nm*x-toDrop->x)/(nm-1);
      y = (nm*y-toDrop->y)/(nm-1);
    }
  // zero it
  toDrop = CountedRef<BaseStar>();
  nm--;
}

const BaseStar *NStarMatch::GetMatch(const unsigned Index) const
{
  if (Index >= stars.size()) return NULL;
  return stars[Index];
}




/************************* NStarMatchList ***********************/


#include <algorithm> // for sort, binary_search
#include <listmatch.h>
#include <starmatch.h>

int NStarMatchList::MatchAnotherList(const BaseStarList &ToMatch, 
				     const double MaxDist, 
				     const BaseStar *EmptyStar,
				     const string &Tag)
{
  emptyStars.push_back(EmptyStar);
  unsigned Index = emptyStars.size() - 1;

  // construct a tag (used in write()) if not provided
  string tag = Tag;
  if (tag == "")
    {
      char chi[80];
      sprintf(chi,"%d",Index+1);
      tag = chi;
    }
  tags.push_back(tag);

  StarMatchList *sml = ListMatchCollect((const BaseStarList &)*this,
					ToMatch,
					MaxDist);
  cout << " matches " << sml->size() << endl;

  vector<const BaseStar *> matchedObjects;

  for (StarMatchIterator  i = sml->begin(); i != sml->end(); ++i)
    {
      StarMatch &sm = *i;
      NStarMatch &nsm = dynamic_cast<NStarMatch &>(*(sm.s1));
      BaseStar &newMatch = *(sm.s2);
      nsm.AddMatch(newMatch,Index);
      matchedObjects.push_back(&newMatch);
    }
  
  /* sort the vector of (pointers to) matches so that we can use a
     binary search when looking for unmatched objects */
  std::sort(matchedObjects.begin(), matchedObjects.end());

  // loop on unmatchedobjects 
  int unmatchCount = 0;
  for (BaseStarCIterator i = ToMatch.begin(); i != ToMatch.end(); ++i)
    {
      const BaseStar *b = *i;
      // was it matched?
      if (binary_search(matchedObjects.begin(), matchedObjects.end(), b))
	continue; /* yes */
      /* no : add it as a single new object */
      NStarMatch *newMatch = new NStarMatch();
      newMatch->AddMatch(*b, Index);
      push_back(newMatch);
      unmatchCount ++;
    }
  cout << " unmatchCount " << unmatchCount << endl;
  int nMatches = sml->size();
  delete sml;
  return nMatches;
}

/* cannot use the default StarList::write, because
   NstarMatch may have "empty matches" (i.e. null pointers), 
   and the way to write them is to substitute emptyStars,
   user provided
*/


void  NStarMatchList::ApplyCountCut(const int MinCount)
{
  for (NStarMatchIterator i = begin(); i != end(); )
    {
      NStarMatch &m = **i;
      if (m.NumberOfMatches() < MinCount)  i = erase(i);
      else ++i;
    }
}

#include <fstream>
#include <iomanip>

int  NStarMatchList::write(const string &FileName) const
{
  ofstream pr(FileName.c_str());
  if (!pr)
    {
      cerr << " StarList::write : could not open " << FileName << endl;
      return 0;
    }
  //  ios::fmtflags  old_flags =  pr.flags(); 
  pr  << resetiosflags(ios::scientific) ;
  // pr  << setiosflags(0) ;
  pr  << setiosflags(ios::fixed) ;
  pr<< setprecision(8);
  string format;
  for (unsigned k=0; k < emptyStars.size(); ++k)
    {      
      format += " "+emptyStars[k]->WriteHeader_(pr,tags[k].c_str());
    }
  /*  dx, dy & dass can be computed from x and y's, but
      when using single precision, the results are poor, so
     "pre" compute them
  */
  // add dx dy dass
  for (unsigned k1=0; k1 < emptyStars.size(); ++k1)
    for (unsigned k2=k1+1; k2 < emptyStars.size(); ++k2)
      {
	pr << "# dx"<< tags[k1] << tags[k2] << " : difference in x" << endl;
	pr << "# dy"<< tags[k1] << tags[k2] << " : difference in y" << endl;
	pr << "# dass"<< tags[k1] << tags[k2] << " : assoc distance" << endl;
      }
  pr << "# format " << format << endl ;
  pr << "# end " << endl ;
  for (NStarMatchCIterator i = begin(); i != end(); ++i)
    {
      const NStarMatch &nsm = **i;
      for (unsigned k=0; k < emptyStars.size(); ++k)
	{
	  const BaseStar *current = nsm.GetMatch(k);
	  if (!current) 
	    current = emptyStars[k];
	  current->writen(pr);
	  pr << " ";
	}
      // add dx dy dass
      for (unsigned k1=0; k1 < emptyStars.size(); ++k1)
	{
	  const BaseStar *c1 = nsm.GetMatch(k1);
	  for (unsigned k2=k1+1; k2 < emptyStars.size(); ++k2)
	    {
	      const BaseStar *c2 = nsm.GetMatch(k2);
	      if (c1 && c2)
		pr << c1->x-c2->x << ' ' << c1->y-c2->y << ' '
		   << c1->Distance(*c2) << ' ';
	      else
		pr << " -1 -1 -1 ";
	    }
	}
      // end of distances
      pr << endl;
    }
  pr.close();
  return 1;
}



#ifdef STORAGE

void NStarMatchList::MatchAnotherList2(const BaseStarList &List, const int rank, const double MaxDist, const BaseStar *EmptyStar)
{
  // here, contrary to MatchAnotherList, we find, for each star of List, 
  // the nearest star in nstarmatchlist and not the opposite.
  // this prevent us to have several matches per star.

  
  emptyStar[rank] = EmptyStar;
  

  
  // create as many lists of stars and associated finders as lists to match
  //(BaseStarList *)*AllStars = new (BaseStarList*)[nlists];
  
  FastFinder **finder;
  finder = new FastFinder*[nlists];
  
  for(int current_rank=0; (current_rank<nlists) && (current_rank!=rank) ; current_rank++) {
    
    // for each rank, create a BaseStarList and a finder for this BaseStarList.
    BaseStarList *AllStars = new BaseStarList();
    for (NStarMatchIterator nsi = begin(); nsi != end(); ++nsi) {
      if (nsi->star[current_rank])
	AllStars->push_back(nsi->star[current_rank]);
    }
    if(AllStars->size()==0) // list is empty
      finder[current_rank] = NULL;
    else
      finder[current_rank] = new FastFinder(*AllStars);
    
  }
  
  // now loop on the input list and do the match
  
  for (BaseStarCIterator bsi = List.begin(); bsi != List.end(); ++bsi) {
    const Point *p1 = *bsi;
    
    double minDist2 = MaxDist*MaxDist*10;
    int best_rank = 0;
    const BaseStar* best_neighbour = NULL;
    
    for(int current_rank=0; (current_rank<nlists) && (current_rank!=rank) ; current_rank++) {
      
      if(!finder[current_rank]) // list is empty
	continue;
      
      Point p2 = this->ApplyTransfo(rank,current_rank,*p1); // p2 is tranformed to the current rank
      const BaseStar* neighbour = finder[current_rank]->FindClosest(p2,MaxDist);
      if(!neighbour) // no match
	continue;
      //cout << "there is a match" << endl;
      // keep the best neighbour comparing with other ranks
      double dist2 = p2.Dist2(*neighbour);
      if(dist2<minDist2) {
	minDist2=dist2;
	best_neighbour=neighbour;
	best_rank=current_rank;
      }
    }
    
    if(best_neighbour) {
      // there is a best match
      // find in which nstarmatch this happened, and save our star
      for (NStarMatchIterator nsi = begin(); (nsi != end()) ; ++nsi) {
	if( nsi->star[best_rank] == best_neighbour) { // it is this one
	  const BaseStar* original = *bsi;
	  nsi->star[rank] = original;
	  nsi->points[rank] = *original;
	  nsi->actualMatches++;
	  break;
	}
      }
    }else{
      // no match so we add this star to the list
      const BaseStar *original = *bsi;
      NStarMatch match(original, rank, nlists);
      push_back(match);
      continue;
    }
  }
  
  // now delete stuff
  for(int current_rank=0; (current_rank<nlists) && (current_rank!=rank) ; current_rank++) {   
    if(finder[current_rank]) {
      if(finder[current_rank]->baselist)
	delete finder[current_rank]->baselist;
      delete finder[current_rank];
      //delete AllStars[current_rank];
    }
  }
  delete[] finder;
  //delete[] AllStars;
}


#endif

