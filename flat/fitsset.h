#ifndef FITSSET__H
#define FITSSET__H

#include <string>
#include <vector>

using namespace std;

/*! \file
   \brief sets of fits files with homogeneous sizes, filter, and chip id

   */




//! container for fits files that have same sizes, filter, and come from the same chip within a mosaic.
class FitsSet
{
  private  :
vector<string> fitsNames;
int nx,ny,nxtot, nytot;
string all_names;

 public :
    //! constructor from a file that contains file names (one per line). 
   /*! Those can be actual filenames or DbImage names. */
 FitsSet(const string &FileName, const bool CheckFilter = true, const bool CheckCCD = true);

#ifndef SWIG
    string operator [] (const unsigned int i) const { return fitsNames[i];}
#endif

    //! Number of file in the set.
    int size() const { return fitsNames.size();}

    //! returns a string containing all file names.
    string AllNames() const { return all_names;}

    //! x size of images (excluding overscan if any).
    int Nx() const { return nx;}

    //! y size of images.
    int Ny() const { return ny;}

    //! total x size of images (including overscan)
    int NxTot() const { return nxtot;}

    //! same for y.
    int NyTot() const { return nytot;}
};





#endif /* FITSSET__H */
