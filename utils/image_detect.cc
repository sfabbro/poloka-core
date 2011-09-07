#include "reducedimage.h"
#include "allreducedimage.h"
#include "detection.h"

static string DetectionsFile(const ReducedImage& Im) { 
  return Im.Dir() + "det.list";
}


struct RunDetection {

  RunDetection() : fixedPositions(false) {}

  BaseStarList Positions;
  bool fixedPositions;
  ReducedImageRef Ref;


  bool operator () (ReducedImageRef Im) const {        

    if (!Im->MakeFits() && !Im->MakeWeight()) {
      cerr << " cannot run detection for " << Im->Name() 
	   << " without both image and weights\n";
      return false;
    }

    string filename = DetectionsFile(*Im);
    DetectionList detections;

    if (Positions.empty() && FileExists(filename)) {
      detections.read(filename);
      return true;
    }
    
    //block to save memory
    {
      DetectionProcess detectionProcess(Im->FitsName(), Im->FitsWeightName(), 
					Im->Seeing(), Im->Seeing());
      
      if (Positions.empty())
	detectionProcess.DoDetection(detections);
      else
	detectionProcess.DetectionScoresFromPositions(Positions, detections, fixedPositions);
    }

    //block to save memory
    {
      /* 
	 fill in the Detection's block that concerns the ref:
	 - flux on the ref at the candidate position (measured in the same conditions
         as the candidate (this is why we use Im->Seeing() rather the Ref->Seeing())
	 - nearest object
      */
      ReducedImageRef ref;
      if (Ref)
	ref = Ref;
      else
	ref = Im;
      DetectionProcess refDet(ref->FitsName(), ref->FitsWeightName(), 
			      Im->Seeing(), Im->Seeing());
      refDet.SetScoresFromRef(detections, *ref);
    }

    detections.write(filename);
    return true;
  }

};


static void usage(const char *progName) {
  cerr << "Usage: " << progName << " [OPTION]...[DBIMAGE]...\n"
       << "Detect and measure point sources using adapted gaussian filter\n\n"
       << "   -r DBIMAGE: specify a reference (for matched det and subs)\n"
       << "   -m FILE   : match all detections into <file>\n"
       << "   -p FILE   : specify a star list for initial positions \n"
       << "   -f        : fixed position measurement (need -p)\n";
}

int main(int nargs, char **args) {

  if (nargs < 2) { usage(args[0]); return EXIT_FAILURE; }

  ReducedImageList imList;
  RunDetection runDetection;
  string matchFile("");

  for (int i=1; i<nargs; ++i) {
    char *arg = args[i];
    if (arg[0] != '-') {
      ReducedImageRef im = ReducedImageNew(arg);
      if (!im || !im->IsValid()) { 
	cerr << " not a valid dbimage: " << arg << endl;
	continue;
      }
      imList.push_back(im);
      continue;
    }

    switch (arg[1]) {
    case 'r': { runDetection.Ref = ReducedImageNew(args[++i]); continue; }
    case 'm': { matchFile = args[++i]; continue; }
    case 'p': { runDetection.Positions.read(args[++i]); continue; }
    case 'f': { runDetection.fixedPositions = true; continue; }
    default: usage(args[0]); return EXIT_FAILURE;
    }
  }

  for_each(imList.begin(), imList.end(), runDetection);

  if (!matchFile.empty() && imList.size() > 1) {
    DetectionList detsOnRef;
    if (runDetection.Ref) {
      runDetection(runDetection.Ref);
      detsOnRef.read(DetectionsFile(*runDetection.Ref));
    } else {
      detsOnRef.read(DetectionsFile(*imList.front()));
    }      
    BaseStarList *positions = Detection2Base(&detsOnRef);
    MatchedDetectionList matchedDetections(detsOnRef);
    for (ReducedImageIterator it=imList.begin(); it != imList.end(); ++it) {
      DetectionList det(DetectionsFile(**it));
      matchedDetections.OneToOneAssoc((*it)->Name(), det);
    }
    matchedDetections.ApplyCuts();
    matchedDetections.write(matchFile);
  }

  return EXIT_SUCCESS;
}
