#ifndef IMAGESUBTRACTION__H
#define IMAGESUBTRACTION__H

#include "reducedimage.h"
#include "psfmatch.h"
#include "fileutils.h"
#include "detection.h" /* for DetectionList */
/*! \file */
 

string SubtractedName(const string &RefName, const string &NewName);
    
//! For subtracting images using the Alard kernel fit technique. A basic assumption: Ref and New are already geometrically aligned.


#include "persistence.h"

class ImageSubtraction : public ReducedImage, public PsfMatch {


  private :

CLASS_VERSION(ImageSubtraction,1);
#define ImageSubtraction__is__persistent

  public :
    //! the constructor takes a copy of both ReducedImage.
    ImageSubtraction(const string &Name, const ReducedImageRef RefImage, const ReducedImageRef NewImage);

    //! as above but uses kernel found by other means. used for fakes SN
    ImageSubtraction(const string &Name, const ReducedImageRef RefImage, const ReducedImageRef NewImage, const PsfMatch *AMatch);

#ifdef STORAGE
    ImageSubtraction(const string &Name);
#endif

    //! Carry out the kernel fit, convolve, and outputs the subtraction image.
    ImageSubtraction(); // necessary for persistence.

    virtual const string  TypeName() const { return "ImageSubtraction";}

    bool MakeFits() ;

    //! 
    bool MakeWeight();

    //! carry out the detection. Default is matched filter plus apreture photometry.
    bool MakeCatalog();
    // ! Dead frame is the OR of Ref and New dead frames.
    bool MakeDead();
    //!  Satur frame is the OR of Ref and New satur frames.
    bool MakeSatur();
    //!  Mask the saturation on the subtraction
    bool MaskSatur();
    //!  Mask null weights on the subtraction
    bool MaskNullWeight();

    bool RunDetection(DetectionList &Detections,
		      const BaseStarList* Positions = NULL);
#ifdef STORAGE
    //!  DeadAndSatur frame is the OR of Ref and New dead and satur frames.
    bool MakeDeadAndSatur();
#endif /* STORAGE */
    //!  Name of the candidates list
    string CandName() const;
    string AllCandName() const;
    string CandCutName() const;
    string CandScanName() const;
    string CandCutScanName() const;

    string DetectionsName() const { return Dir()+"det.list";}
    string MatchedCatName() const { return Dir()+"matcheddet.list";}
    
    //!
    ReducedImage *Clone() const;

    string AllCandidateCatalogName() const { return AddSlash(Dir())+"allcand.list";}

    ~ImageSubtraction();

#ifndef SWIG
    ClassDef(ImageSubtraction,1);
#endif /* SWIG */
};

typedef CountedRef<ImageSubtraction> ImageSubtractionRef;

#endif /* IMAGESUBTRACTION__H */
