#include <iomanip>
#include "fileutils.h"
#include "detection.h"
#include "nstarmatch.h"

// uncomment CHECK_MULTI_MATCH to check multimatches
// that cause pbs when swaping ra,dec/x,y
//  pointers to stars are saved in the list
// so swaping happens several times when star is 
// several times in the list.
// this is supposed to be fixed now in nstarmatch
// so this tested is not needed any more
// (it slows down computation)

//#define CHECK_MULTI_MATCH

 void usage(const char *progName)
{
  cerr << progName << " [-o <merged_list> -d <x> -s <x>]  <candlist(s)> \n"
       << "  match candidates list (in det.list format) \n"
       << "   -o <merged_list> output for the matched list (default ./matchradec.list)\n"
       << "   -d <x> max association distance in arcsec (default to 1) \n"
       << "   -s <x> min combined signal to noise (default to 3) \n";
  exit(-1);
}

 void swap_xy_radec(Detection *star)
{
  if(star==NULL)
    return;
  double xtmp = star->x;
  double ytmp = star->y;
  star->x = star->Ra();
  star->y = star->Dec();
  star->Ra() = xtmp;
  star->Dec() = ytmp;
  // should do the same for variance and all that, but useless for now
}

 void RaDecMatch(const vector<string> &ToMatch, const string &MatchedListName, 
		       const double &MaxDist, const double &MinStn)
{
  int n = ToMatch.size();
  NStarMatchList mList(n);

  vector <DetectionList*> detList;//to delete aftewards.
  GtransfoIdentity ident;
  
  for (int i=0; i<n ; ++i) {
    DetectionList *current_list = new DetectionList(ToMatch[i]);
    
    detList.push_back(current_list);
    cout << " Matching " << ToMatch[i] << " (" << current_list->size() << " stars) ... \n";
    
    // swap between x y and ra dec to match in ra dec using fastfinder
    
    
    for (DetectionIterator di = current_list->begin(); di != current_list->end(); ++di)
      {
	swap_xy_radec(*di);
      }
    
    // Images are already aligned in ra and dec so transfo is identity
    
    mList.SetTransfo(&ident, i, 0);

    mList.SetTransfo(&ident, 0, i);
    
    //mList.MatchAnotherList(*Detection2Base(current_list), i, MaxDist, new Detection());
    mList.MatchAnotherList2(*Detection2Base(current_list), i, MaxDist, new Detection());
    //cout << " done" << endl;
  }
  
  // now re-swap xy ra dec. Be careful because MatchAnotherList set everything 
  // to zero on the non-matched stars
  int deleted = 0;
  
#ifdef CHECK_MULTI_MATCH
  vector< Detection *> swap_done;
  cout << "Loop on list and swap (ra,dec)/(x,y) ..." << endl;
  int count=0;
#endif
  int ncount = mList.size();

  for (NStarMatchIterator mi = mList.begin(); mi != mList.end(); )
    {
#ifdef CHECK_MULTI_MATCH
      count ++;
      if(count%100==0)
	cout << count << "/" << ncount << endl;
#endif
      double stn = 0;
      double x = 0;
      double y = 0;
      for (int i=0; i< mi->npointers; ++i)
	{
	  if (mi->StarExist(i))
	    {
	      Detection *current = (Detection*) mi->star[i];
#ifdef CHECK_MULTI_MATCH
	      bool isdone=false;
	      for(unsigned int is=0;is<swap_done.size();is++) {
		if(swap_done[is]==current) {
		  cout << "is done already = this is due to a bug in nstarmatch" << endl;
		  isdone = true;
		  break;
		}
	      } 
	      if(isdone)
		continue;
	      swap_done.push_back(current);
#endif
	      swap_xy_radec(current);
	      
	      if (x == 0)
		{
		  x = current->x;
		  y = current->y;
		}
	    }
	}
      
      
      for (int i=0; i< mi->npointers; ++i)
	{
	  if (mi->StarExist(i)) {
	    Detection *current = (Detection*) mi->star[i];
	    stn += current->Sig2Noise()*current->Sig2Noise();
	  }else{	    
	    mi->star[i] = new Detection(x,y,0);	      
	  }
	}
      if (sqrt(stn) < MinStn) 
	{
	  mi = mList.erase(mi);
	  deleted++;
	}
      else 
	{
	  ++mi;
	}
    }
  
  cout << " Final matched list size : " << mList.size() << "/" << ncount << " stars " << endl;
  cout << "Deleted "  << deleted << endl;
  // write the list with better precision than usual, just for ra and dec
  ofstream mout(MatchedListName.c_str());
  mout << setiosflags(ios::fixed);
  mout << setprecision(5);
  mList.write(mout);

  for (int i=0; i<n ; ++i) delete detList[i];
}

int main(int nargs, char **args)
{

  if (nargs < 3){usage(args[0]);}

  // default values
  string matchedListName = "matchradec.list";
  vector<string> toMatch;
  double maxdist = 1;
  double minstn = 3;

  // decode arguments
  for (int i=1; i<nargs; ++i)
    {
      char *arg = args[i];
      if (arg[0] != '-') 
	{
	  if (!FileExists(arg))
	    {
	      cerr << " No detection list " << arg << endl;
	      continue;
	    }
	  toMatch.push_back(arg);
	  continue;
      }
      switch(arg[1])
	{
	case 'o' : ++i; matchedListName = args[i]; break;
	case 'd' : ++i; maxdist = atof(args[i]); break;
	case 's' : ++i; minstn = atof(args[i]); break;
	default : usage(args[0]);
	}
    }  

  // do the match
  RaDecMatch(toMatch, matchedListName, maxdist/3600., minstn);
  
  return EXIT_SUCCESS;
}
