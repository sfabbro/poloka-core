#include <iostream>
#include <fstream>

#include "gtransfo.h"
#include "nstarmatch.h"
#include "basestar.h"
#include "vutils.h" /* for DArrayMedian */
#include "fastfinder.h"


NStarMatch::NStarMatch(int n)
{
  npointers=n; actualMatches = 1;
  points= new Point[npointers];
  star = new const BaseStar*[npointers];
  for (int i=0; i<n; i++) star[i]=NULL;
}

NStarMatch::~NStarMatch()
{
  delete [] points;
  delete [] star;
}

NStarMatch::NStarMatch(const NStarMatch &O) : BaseStar(*this)
{
  npointers = O.npointers;
  actualMatches = O.actualMatches;
  points= new Point[npointers];
  memcpy(points, O.points, npointers*sizeof(Point));
  star = new const BaseStar*[npointers];
  memcpy(star, O.star, npointers*sizeof(BaseStar*));
}

NStarMatch::NStarMatch(const BaseStar *Star, const int i, const int n) : BaseStar(*Star)
{
  npointers = n;
  actualMatches = 1;
  //ref=0;
  points= new Point[n];
  star = new const BaseStar*[n];
  for (int k=0; k<n; k++) star[k]=NULL;
  points[i]=(*Star);
  star[i]=Star;
}

/************************* NStarMatchList ***********************/

NStarMatchList::NStarMatchList(const int n)   : nlists(n),nused(n,n),order(n,n),chi2(n,n),transfo(n,n,NULL)
{
  emptyStar = new const BaseStar*[nlists]; 
  for (int i=0; i< nlists; ++i) emptyStar[i] = NULL;
}

						   

Point NStarMatchList::ApplyTransfo(const int i, const int j, Point pt)
{
  if (!(this->transfo(i,j)))
    {
      //      Gtransfo transformation=GtransfoCompose(*this->transfo(j,0),*this->transfo(0,i));
      SetTransfo(GtransfoCompose(this->transfo(j,0),this->transfo(0,i)),i,j);
    }
    return this->transfo(i,j)->apply(pt);
}

void NStarMatchList::MatchAnotherList(const BaseStarList &List, const int rank, const double MaxDist, const BaseStar *EmptyStar)
{
  FastFinder finder(List);
  emptyStar[rank] = EmptyStar;
  // this is complicated :
  // We want to identify objects in List which will remain unmatched.
  // since List is const, we store matchedObjects in a copy and find unmatched   
  // objects at the end. We cannot work on a non const copy of List 
  // since star pointers of NStarMatch's have to point to List, 
  // and not to a copy.
  // this is better :
  // we just save the values of const BaseStar pointers
 
  //BaseStarList matchedObjects;
  vector<const BaseStar*> matchedObjects;
  for (NStarMatchIterator si = this->begin(); si != this->end(); ++si)
    {
      for (int i =0; (i<this->nlists) && (i!=rank) ; i++)
	{
	  if (!(si->star[i])) continue;
	  const Point *p1 = si->star[i];
	  Point p2 = this->ApplyTransfo(i,rank,*p1);
	  const BaseStar *neighbour = finder.FindClosest(p2,MaxDist);
	  if (!neighbour) continue;
	  bool ismatched = false;
	  for (unsigned int ib=0;ib<matchedObjects.size();ib++) {
	    if(neighbour==matchedObjects[ib]) {
	      ismatched = true;
	      break;
	    }
	  }
	  if(ismatched) {
	    // uncomment this and check if you have multi matches.
	    // in any case, you'd better use the next routine  MatchAnotherList2
	    //cout << "whao it was already matched !!" << endl;
	    continue;
	  }

	  // store the match to identify later leftovers.
          //matchedObjects.push_back(new BaseStar(*neighbour));
          matchedObjects.push_back(neighbour);
	  
	  // handle the match
	  si->star[rank] = neighbour;
	  si->points[rank] = *neighbour;
	  si->actualMatches++;
	  break;
	}
    }
  
//   FastFinder matchedFinder(matchedObjects);
//   for (BaseStarCIterator si = List.begin(); si != List.end(); ++si)
//     {
//       const BaseStar *original = *si;
//       const BaseStar *copy = matchedFinder.FindClosest(*original, MaxDist/100.);
//       if (!copy) // original was never matched
// 	{
// 	  NStarMatch match(original, rank, nlists);
// 	  push_back(match);
// 	}
//     } 
  
  // now we compare values of pointers to check if a detection object
  // is already matched
  for (BaseStarCIterator si = List.begin(); si != List.end(); ++si)
    {
      bool ismatched = false;
      const BaseStar *original = *si;
      for (unsigned int ib=0;ib<matchedObjects.size();ib++) {
	if(original==matchedObjects[ib]) {
	  ismatched = true;
	  break;
	}
      }
      if (!ismatched) // original was never matched
	{
	  NStarMatch match(original, rank, nlists);
	  push_back(match);
	}
    }
}


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


void NStarMatchList::short_write(ostream & pr) const
{
  NStarMatchCIterator itu= begin();
  NStarMatch Nstarm = *itu ; 
  GtransfoIdentity identity;
  if ( itu == end() )
    {
      cerr << " Pas d'Ecriture de StarMatchList vide " << endl ;
      return ;
    }

  for (int i =0; i<this->nlists ; i++)
    {
      char chi[80];
      sprintf(chi,"%d",i+1);
      pr<< "# x"<< chi << " : x of star "<< chi <<endl;
      pr<< "# y"<< chi << " : y of star "<< chi <<endl;
      pr<< "# f"<< chi << " : flux of star "<< chi <<endl;
    }
  pr << "# end " << endl ;
 
  GtransfoIdentity id ;
  for (NStarMatchCIterator it= begin(); it!= end(); it++ )
    {
      const NStarMatch &Nstarm = *it ;

      if (!(Nstarm.StarExist(0)))
	pr << "0.0 0.0 0.0";
      else
	pr << Nstarm.x(0)<< " " << Nstarm.y(0)<< " " << Nstarm.flux(0);
      for (int i =1; i<this->nlists ; i++)
	{
	  if (!(Nstarm.StarExist(i)))
	    pr << " 0.0 0.0 0.0";
	  else
	    pr << " " << Nstarm.x(i)<< " " << Nstarm.y(i)<< " " << Nstarm.flux(i);
	}
      pr << endl ;
    }

}

void NStarMatchList::short_write(const string &filename) const
{
  ofstream pr(filename.c_str()) ;
  short_write(pr);
  pr.close();
}

void NStarMatchList::write(ostream & pr) const
{
  NStarMatchCIterator itu= begin();
  if ( itu == end() )
    {
      cerr << " Cannot write an empty NStarMatchList " << endl ;
      return ;
    }

  for (int i=0; i<nlists; ++i)
    {
      char chi[80];
      sprintf(chi,"%d",i+1);
      emptyStar[i]->WriteHeader_(pr,chi);
    }

   for (int i =0; i<this->nlists ; i++)
     for (int j = i+1; j<this->nlists ; j++)
	 pr << "# dass"<< i << j <<" : association distance for "
	    <<i<<" and "<<j<< endl ; 

   // 
   pr << "# end " << endl ;
 
   int count = 0;
   for (NStarMatchCIterator it= begin(); it!= end(); it++ )
     {
       const NStarMatch &Nstarm = *it ;

       for (int i =0; i<this->nlists ; i++)
	 {
	   const BaseStar *current = Nstarm.star[i]; 
	   if (!current) current = emptyStar[i];
	   current->writen(pr);
	   pr << " " ;
	 }
       for (int i =0; i<this->nlists ; i++)
	 for (int j = i+1; j<this->nlists ; j++) {
	   if(Transfo(i,j))
	     pr << Nstarm.Distance(i,j, *Transfo(i,j))<< " " ;
	   else
	     pr << 0 << " " ;
	 }
      pr << endl ;
      ++count;
    }

}

void NStarMatchList::write(const string &filename) const
{
  ofstream pr(filename.c_str()) ;
  write(pr);
  pr.close();
}

