#include "imagelist.h"
#include "iostream" // for debugging

#ifdef USE_ROOT

#include <TObject.h>
#include <TString.h>

template <class T> void ImageList<T>::Streamer(TBuffer &R__b)
{
  
   if (R__b.IsReading()) {
      ImageList<T>::Class()->ReadBuffer(R__b, this);
      {
	cout << " reading an image list " << endl;
         clear();
         int R__i, R__n;
         R__b >> R__n;
         for (R__i = 0; R__i < R__n; R__i++) {
	   TString R__str; 
	   R__str.Streamer(R__b); 
	   string imageName = R__str.Data(); 
	   push_back(ReducedImageRead(imageName));
         }
      }
   } else {
      ImageList<T>::Class()->WriteBuffer(R__b, this);
      {
	cout << " writing an image list " << endl;
         R__b << int(size());
         ImageList<T>::iterator R__k;
         for (R__k = begin(); R__k != end(); ++R__k) {
	   ReducedImage *i = *R__k;
           TString R__str = i->Name().c_str();
	   R__str.Streamer(R__b);
         }
      }
   }
}


#endif /*USE_ROOT */
