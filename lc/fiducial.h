// This may look like C code, but it is really -*- C++ -*-
#ifndef FIDUCIAL__H
#define FIDUCIAL__H

#include <iomanip>

#include <countedref.h>
#include <reducedimage.h>
#include <iomanip>

//! a template to use when an pointer element belongs to an image
template<class S> 
class Fiducial : public S {

protected:

  CountedRef<ReducedImage> rim;
  
public:

  Fiducial() {}

  Fiducial(const S *Fid) : S(*Fid) {}

  Fiducial(const ReducedImage *Rim) : rim(Rim) {}

  Fiducial(const S *Fid, const ReducedImage *Rim) : S(*Fid), rim(Rim) {}

  const ReducedImage* Image() const { return rim; }

  bool HasImage(const ReducedImage *Rim) const
  {
    return (rim->Name() == Rim->Name());
  }

  void AssignImage(const ReducedImage *Rim) 
  { 
    if (rim) { cerr << " Fiducial::AssignImage() : Error : " 
		    << rim->Name() << " already assigned \n";  return; }
    rim = Rim;
    //cout << "in AssignImage rim=" << rim->Name() << endl;
  }

  ostream& write_header(ostream &Stream=cout) const
  {
    if (rim) 
      {
	Stream << "# name : DbImage name \n";
	Stream << "# band : integer for the band (ascii code) \n";
      }
    S::write_header(Stream);
    Stream << "# end \n";

    return Stream;
  }

  friend ostream& operator << (ostream &Stream, const Fiducial<S>& Fs)
  {
    ios::fmtflags oldflags = Stream.flags();
    Stream.setf(ios::fixed);
    if (Fs.rim)
      Stream << setw(12) << setprecision(2) << Fs.rim->JulianDate()
	     << int(Fs.rim->Band()[0]);
    Stream << (S) (Fs);
    Stream.flags(oldflags);    
    return Stream;
  }

  friend bool IncreasingJulianDate(const Fiducial<S> *one, const Fiducial<S> *two)
   { return (one->rim->JulianDate() < two->rim->JulianDate()); }

  friend bool IncreasingSeeing(const Fiducial<S> *one, const Fiducial<S> *two)
   { return (one->rim->Seeing() < two->rim->Seeing()); }

};


#endif // FIDUCIAL__H
