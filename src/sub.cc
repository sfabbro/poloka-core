#include <iostream>
#include <cmath>

#include "sub.h"
#include "fitsimage.h"
#include "wcsutils.h"
#include "imagematch.h"
#include "gtransfo.h"
#include "subimage.h"
#include "imagesum.h"
#include "imagesubtraction.h"
#include "reducedutils.h"

/*! \page newsubfile Syntax of the "subfile"
a subfile is a very simple text file that describes what should be 
subtracted from what. We define 3 set of images, labelled "ref", "new1" and "new2",
where "new1 and "new2" refer to the same epoch. We can either run a

 - "split subtraction" in which case, we compute new1-ref, new2-ref and
new1+new2-ref; the detection will request to find an object in
the 3 subtractions (with different significances). 
 - or a single subtraction in which case new1 and new2 groups are merged,
and we compute new-ref.

  Before subtracting each group of images is aligned and coadded (using
ImagesAlignAndSum). At the moment, seeings in a given stack are 
NOT equalized, because the subtraction
makes use of a kernel fitting (see KernelFit) that can handle "mixed PSF's".

A "subfile" details the contents of these subtraction components through
DbImage names. The choice between a single subtraction and "split" subtraction
is specified by putting "ONE SUBTRACTION" in this file. "ADDFAKES" is a 
statement that drives the code to add fake supernovae in the input (new) images
(they are not altered anyway!) for detection efficiency studies. This option
roughly doubles the processing time. Blank lines in the subfile are ignored.
By default, newmake_sub reads a file named subfile in the current directory.

"FIXGEO" aims at fixing the geometric reference used to align all the images.
Here is an example of a "subfile":

\code
# subtraction for fieldP1ccd9
REF
502809o09 <SubImage>
502810o09
502811o09
502813o09
502814o09
502815o09
502816o09
502817o09
NEW1
504604o09
504605o09
504606o09
504607o09
NEW2
504613o09
504614o09
504615o09
504616o09
#ADDFAKES
ONE SUBTRACTION
#FIXGEO
502809o09
\endcode

In this specific example, the split of images between NEW1 and NEW2 is
irrelevant since we run a "simple" subtraction (as specified by ONE SUBTRACTION). ADDFAKES is commented and no fakes will be added.


*/

typedef enum TypeSub{TypeRef =0, TypeNew1, TypeNew2, TypeNew};

struct DatSim {
  int numberOfFakes;
  double minMag, maxMag;
  
  DatSim() { numberOfFakes = 100; minMag = 22; maxMag = 26;}
  void LitDataCards(DataCards &);
  DatSim(const string &FileName);
  void Print();
};


#ifndef M_PI
#define  M_PI           3.14159265358979323846  /* pi */
#endif

#include "stdlib.h"
// does not rename: builds an alternate name
static void RenameDbImage(const DbImage& D, const string &NewName)
{
  char command[256];
  sprintf(command,"cd %s/..; ln -s %s %s\n",D.Dir().c_str(), 
	  D.Name().c_str(), NewName.c_str());
  if (system(command)!=0)
    {
      cerr << " command : " << endl
	   << command 
	   << " failed " << endl;
    }
}


static const char* last_word(const char *line, const char sep = ' ')
{
  const char *p = line + strlen(line) - 1;
  // skip trailing sep's
  while ( (*p == sep) && (p > line)) p--;
  while ( (*p != sep) && (p > line)) p--;
  return p;
}

static string first_word(const char *line, const char sep = ' ')
{
  int start;
  int end_line = strlen(line);

  while ( line[start] == sep && start < end_line) start++;
  int end = start;
  while ( line[end] != sep &&  end < end_line) end++;
  return string(line).substr(start,end);;
}



static double sqr(double x){return(x*x);}

// Read the subfile. If you change the syntax here, update the example just above.
Sub::Sub(const string &FileName, const bool Overwrite, const bool OnlyDet) : overwrite(Overwrite)
{
  // Boolean for the simulation
  //  FakeList = NULL;

  AssociateGal = false;
  onlyDet = OnlyDet;
  detectOnAllSub = false;
  AddFakes = false;
  FixRef = false;
  StringList ToExtract;
  ImageNameToExtract = "";
  onlyOneSub = false;

  
  FILE *file = fopen(FileName.c_str(),"r");
  if (!file)
    {
      cerr << " ERROR : Sub::Sub : cannot open \"" << FileName << "\"" << endl;
      exit(1);
      return;
    }
  char line[512];
  bool inRef = false;
  bool inNew = false;
  bool inFixRef = false;
  NewStack *currentNewStack = NULL;


  while (fgets(line,512,file))
    {
      if (line[0] == '#') continue;
      if (line[strlen(line)-1] == '\n') line[strlen(line)-1] = '\0';
      if (strstr(line,"REF"))
	{ inRef = true; inNew = false; inFixRef = false;
	continue;}
      if (strstr(line,"NEW")) 
	{ 
	  inNew = true; inRef = false;  inFixRef = false;
	  AllNew.push_back(NewStack());
	  currentNewStack = &AllNew.back();
	  string currentNewName = StringToLower(last_word(line)); 
	  // remove trailing white space
	  RemovePattern(currentNewName," ");
	  currentNewStack->name = currentNewName;
	  continue;
	}
      if (strstr(line,"ONE SUBTRACTION"))
	{ onlyOneSub = true; continue;}
      if (strstr(line,"DETECT ON ALL SUB"))
	{ detectOnAllSub = true; continue;}
      if (strstr(line,"ADDFAKES")) 
	{ AddFakes = true; cout << "ADDFAKES " << AddFakes << endl; continue;}

      if (strstr(line,"ASSOCIATEGAL"))
	{ AssociateGal = true; continue;}
      if (strstr(line,"FIXGEO"))// dont put REF in that name!!
	{
	  inNew = false ; inRef = false; inFixRef = true;
	  FixRef= true;
	  continue;
	}

      if (!inRef && !inNew && !inFixRef)
	{
	  cerr << " ERROR : unrecognised syntax in file " << FileName << endl ;
	  cerr << line << endl ;
	  cerr <<" stop here " << endl;
          exit(1);
	  return;
	}
      char *start_line = line+strspn(line," "); // skip spaces
      if (start_line[strlen(start_line)-1] == '\n' )
	{
	  start_line[strlen(start_line)-1] = '\0' ;
	}
      if (strlen(start_line) == 0) continue; // skip blank lines
  
      string currentName = first_word(start_line);
      RemovePattern(currentName," ");
      ReducedImage current(currentName);
      if (!current.IsValid())
	{
	  cerr << " cannot find DbImage : \"" << currentName << "\"" << endl;
          if (onlyDet) 
	    { 
	      cerr << " tolerable since you asked for detection only" << endl;
	      continue;
	    }
	  exit(1);
	}

      if (inFixRef)
	{
	  GeomRefName = currentName;
	  cout << " The Geometric reference is fixed : " 
	       << GeomRefName << endl; 

	}
      if (strstr(line, "SubImage") != 0)
	{
	  ToExtract.push_back(currentName);
	}
      if (inRef)
	{
	  Ref.push_back(currentName);
	}
      else if (inNew)
	{
	  currentNewStack->push_back(currentName); 
	}      
      AllInputImages.push_back(currentName);
      //cerr << line << endl;
    }
  if (AllInputImages.size() == 0)
    {
      cerr << " ERROR : Did not find any image name in " <<  FileName << " : stop here " << endl;
      exit(1);
    }
  fclose(file);
  if (ToExtract.size() > 1)
    {
      cerr << " ERROR : We cannot extract subimages from more than 1 image " << endl
	   << " and get aligned images " << endl;
      cerr << " you have requested to extract subimages from : " << endl
	   << ToExtract << endl;
      exit(1);
    }
  if (ToExtract.size() > 0 )
    {
      ImageNameToExtract = ToExtract.front();
      if (FixRef && ImageNameToExtract != GeomRefName)
	{
	  cerr << " ERROR : We cannot extract a sub image from " 
	       << ImageNameToExtract 
	       << endl 
	       << "and impose a different geometric ref " 
	       << GeomRefName << endl;
	  exit(1);
	}
    }
  if (!CheckNewNames())
    {
      cerr << " ERROR : Giving up " << endl;
      exit(1);
    }
  cout << " Read subfile successfully " << endl;
}

int Sub::CheckNewNames()
{
  if (AllNew.size() == 1) return 1;
  for (unsigned i=0; i< AllNew.size(); ++i)
    {
      if (AllNew[i].Name() == "new")
	{
	  cerr << " \"new\" is a reserved name for the global new image" 
	       << endl;
	  return 0;
	}
        for (unsigned j=i+1; j< AllNew.size(); ++j)
	  {
	    if (AllNew[i].Name() == AllNew[j].Name())
	      {
		cerr << " you gave identical names to new images : " 
		     << '\"' << AllNew[i].Name() << "\" and \"" 
		     << AllNew[j].Name() << "\"" << endl;
		cerr << " relabeling new images as new%d " << endl;
		for (unsigned k = 0; k < AllNew.size(); ++k)
		  {
		    char name[16];
		    sprintf(name,"new%d",k+1);
		    AllNew[k].name = string(name);
		  }
		return 1;
	      }
	  }
    }
  return 1;
}
	  


void Sub::RemoveImage(const string &BannedName)
{
  for (unsigned i=0; i< AllNew.size(); ++i) AllNew[i].Remove(BannedName);
  AllInputImages.Remove(BannedName);
  Ref.Remove(BannedName);
}

StringList Sub::AllNewNames() const
{
  StringList res;
  for (unsigned i=0; i< AllNew.size(); ++i)
    {
      const StringList &current = AllNew[i];
      for (StringCIterator k = current.begin(); k != current.end(); ++k)
	res.push_back(*k);
    }
  return res;
}
	     

void Sub::SubstituteName(const string &Original, 
			 const string &Substitution)
{
  for (unsigned i=0; i< AllNew.size(); ++i) 
    AllNew[i].Substitute(Original, Substitution);
  AllInputImages.Substitute(Original, Substitution);
  Ref.Substitute(Original, Substitution);
}


ReducedImage* Sub::ExtractSubimage(const string &SubName)
{
  /* first compute the union of Frame(largeImage)*Frame(anyOtherImage) */
  ReducedImage large(ImageNameToExtract);
  FitsHeader largeHead(large.FitsName());
  Frame largeFrame(largeHead, WholeSizeFrame);
  Gtransfo *largePix2RaDec;
  if (!WCSFromHeader(largeHead, largePix2RaDec))
    {
      cerr << " ERROR : cannot handle a large reference without a WCS " 
	   << endl;
      exit(1);
    }
  Gtransfo *largeRaDec2Pix = 
    largePix2RaDec->InverseTransfo(0.1 /* accuracy in pixels*/, largeFrame);
  delete largePix2RaDec;

  Frame toExtract; // =(0,0,0,0)
  for (StringIterator i = AllInputImages.begin() ; i != AllInputImages.end(); )
    {
      const string &currentName = *i;
      if (currentName == ImageNameToExtract) {++i;continue;}
      ReducedImage current(currentName);
      FitsHeader currentHead(current.FitsName());
      CountedRef<Gtransfo> current2Large;
      if (HasLinWCS(currentHead))
	{
	  Gtransfo *currentPix2RaDec;
	  WCSFromHeader(currentHead, currentPix2RaDec);
	  current2Large = GtransfoCompose(largeRaDec2Pix, currentPix2RaDec);
	  delete currentPix2RaDec;
	}
      else // go the hard way
	{
	  CountedRef<Gtransfo> large2Current;
	  if (ImageListMatch(current, large, current2Large, large2Current))
	    {
	      cerr << " could not match " << current.Name() << " and " << large.Name() << endl;
	      cerr << " we forget " <<  current.Name() << endl;
	      i = AllInputImages.erase(i);
	      RemoveImage(current.Name());	      
	      continue;
	    }
	}
      cout << " transfo between " << large.Name() << ' ' << " and " 
	   << current.Name() << endl 
	   << *current2Large;
      Frame currentFrame(currentHead, WholeSizeFrame); // 
      currentFrame = currentFrame.ApplyTransfo(*current2Large);
      /* if we are at the first image in the loop we have to initialize
	 "toExtract", else we just "increment" it so that we have the union
	 of all input images at the end */
      if (toExtract.Area() == 0.)
	toExtract = (largeFrame*currentFrame);
      else
	toExtract += (largeFrame*currentFrame);
      cout << " after adding " << current.Name() 
	   << ", extract frame = " << toExtract;
      ++i;
    }
  delete largeRaDec2Pix;
  SubImage *subImage = new SubImage(SubName, ImageNameToExtract, toExtract);
  // to be able to match, "make" image and catalog
  subImage->MakeFits(); // need image to make catalog!
  subImage->MakeWeight(); // need also weight to make catalog!
  subImage->MakeCatalog();
  return subImage;
}



void Sub::FindGeometricReference()
{  
  if (ImageNameToExtract != "") // the subimage is gone be the geomref
    {
      string subName;
      if (Ref.size() == 1 && Ref.Contains(ImageNameToExtract)) subName = "ref";
      else subName = "sub_"+ImageNameToExtract;
      GeometricReference = ExtractSubimage(subName);
      //  cout << " Ref before subs " << Ref << endl;
      SubstituteName(ImageNameToExtract, subName);
      // cout << " Ref after subs " << Ref << endl;
    }
  else if (GeomRefName != "") // easy
    {
      GeometricReference = new ReducedImage(GeomRefName);
    }
  else
    {
      /* locate the worst seeing image among all */
      ReducedImageList allImages;
      for (StringCIterator i= AllInputImages.begin(); 
	   i!= AllInputImages.end(); ++i)
	allImages.push_back(new ReducedImage(*i));
      
      GeometricReference = allImages.front();
  
      for (ReducedImageCIterator ri = allImages.begin(); 
	   ri != allImages.end(); ++ri)
	{  
	  const ReducedImage *current = *ri;
	  //cout << current->Name() << endl;
	  if (current->Seeing()*current->PixelSize() 
	      > GeometricReference->Seeing()*GeometricReference->PixelSize()) 
	    GeometricReference = current;
	}
    }
  cout << " Choosing " << GeometricReference->Name() << " as GeometricReference " << endl;
}


ReducedImage* Sub::DoOneStack(StringList Components, const string &StackName, int ToDo)
{
  // a single image stack to do with stackname = single image name
  if (Components.size() == 1 && Components.front() == StackName)
    {
      if (StackName != GeometricReference->Name())
	{
	  cerr << " Error : we are trying to make a stack named \"" 
	       << StackName << "\"" << endl
	       << " from " << Components << endl
	       << " And this is not the GeomRef : cannot align with in==out" 
	       << endl;
	  return NULL;
	}
      GeometricReference->Execute(ToDo);
      return GeometricReference;
    }
  // Do the stack
  ReducedImageList unalignedImages(Components);
  ImageSum *stack = ImagesAlignAndSum(unalignedImages, *GeometricReference,
				      StackName, ToDo, RefStack);
  return stack;
}


int Sub::DoIt()
{
  FindGeometricReference();
  ReducedImageList NewStacks;
  
  //  stack the subparts:
  //int toDo = DoFits | DoCatalog | DoDead | DoSatur;
 
  // The dead is not used anymore (the dead are clipped)
  int toDo = DoFits | DoCatalog | DoSatur | DoWeight;
  cerr << "Processing ref stack " << endl ;
  RefStack = DoOneStack(Ref, "ref", toDo);
  if (!RefStack)
    {
      cerr << " ERROR: No RefStack, stop here " << endl;
      return 0;
    }
  for (unsigned i=0; i<AllNew.size(); ++i)
    {
      NewStack &stack = AllNew[i];
      cerr << "Processing " << stack.Name() << endl ;
      stack.newStack = DoOneStack(stack, stack.Name() , toDo);
      // build the subtraction name
      string newName = stack.Name();
      RemovePattern(newName, "new");
      string subName = "sub"+newName;
      stack.sub = new ImageSubtraction(subName, *RefStack, *stack.newStack);
      stack.sub->Execute(DoFits+DoWeight);
      if (detectOnAllSub) stack.sub->MakeCatalog();
    }
  // several new stacks
  if (AllNew.size() > 1)
    {
      ReducedImageList allAlignedNew(false);
      for (unsigned i=0; i< AllNew.size(); ++i)
	{
	  ReducedImage *stack = AllNew[i].newStack;
	  ImageSum *sum = dynamic_cast<ImageSum *>(stack);
	  if (sum) // this is an actual sum
	    {
	      ReducedImageList sumComponents = sum->Components();
	      for (ReducedImageIterator c = sumComponents.begin(); 
		   c != sumComponents.end(); ++c)
		allAlignedNew.push_back(*c);
	    }
	  else // it should a single image that is the geomref
	    {
	      allAlignedNew.push_back(stack);
	    }
	}
      GlobalNew = new ImageSum("new",allAlignedNew, RefStack);
      GlobalNew->Execute(DoFits + DoCatalog);
      GlobalSub = new ImageSubtraction("sub", *RefStack, *GlobalNew);
      GlobalSub->Execute(DoFits+DoWeight);
    }
  else // only one subtraction
    {
      // we only have one new/sub
      GlobalSub = AllNew[0].sub;
      if (GlobalSub->Name() != "sub")
	RenameDbImage(*GlobalSub,"sub");
      if (AllNew[0].newStack->Name() != "new")
	RenameDbImage(*AllNew[0].newStack, "new");
    }
  RunDetection();
  return 1;
}

int Sub::ExpectedMagLim()
{
  double nsigma = 3.0;
  int Nimages = AllNew.size() + 1;
  static double* photomRatio = new double[Nimages];
  static double* sigmaSeeing = new double[Nimages];
  static double* sigmaSky = new double[Nimages];
  double ZPT;

  ReducedImage Ref(Ref.front());
  ReducedImageList unalignedNewImages(&AllNew);

  photomRatio[0] = 1.0;
  sigmaSeeing[0] = Ref.Seeing();
  sigmaSky[0] = Ref.SigmaBack();
  ZPT = Ref.AnyZeroPoint();

  int j=1;
  for (ReducedImageIterator i = unalignedNewImages.begin(); i != unalignedNewImages.end(); ++i)
    {
      ReducedImage *ri = *i;
      double err;      
      photomRatio[j] = QuickPhotomRatio(*ri,Ref,err);
      sigmaSeeing[j] = ri->Seeing();
      sigmaSky[j] = ri->SigmaBack();
      j++;
    }

  double sigmaTot2 = 0;
  for (unsigned i=0; i<Nimages; ++i)
    {
      sigmaTot2 += 4*M_PI*photomRatio[i]*sigmaSeeing[i]*sigmaSeeing[i]*sigmaSky[i]*sigmaSky[i];
    }

  cout <<" Limiting Mag (" << nsigma <<", ZeroPoint = " << ZPT << " ) = " << ZPT - 2.5*log(nsigma*sqrt(sigmaTot2)) << endl;

  return 1;
}


void Sub::RunDetection()
{
  DetectionList detsOnGlobal;
  GlobalSub->RunDetection(detsOnGlobal);
  if (AllNew.size() == 1) return; // no need to redo the same thing
  BaseStarList *positions = Detection2Base(&detsOnGlobal);
  MatchedDetectionList matchedDetections(detsOnGlobal);
  for (unsigned k=0; k<AllNew.size(); ++k)
    {
      ImageSubtraction* partialSub = AllNew[k].sub;
      DetectionList partialDet;
      partialSub->RunDetection(partialDet, positions);
      matchedDetections.OneToOneAssoc(partialSub->Name(), partialDet);
    }
  matchedDetections.ApplyCuts();
  matchedDetections.write(GlobalSub->MatchedCatName());
}
      



Sub::~Sub()
{
  /*  if (FakeList) delete FakeList; */
}


NewStack::~NewStack()
{
}


#ifdef STORAGE



static void ReducedToTransformed(const ReducedImageList &RIList, ReducedImageList &TIList)
{
  TIList.clear();
  for (ReducedImageCIterator rii = RIList.begin(); rii != RIList.end(); ++rii)
    {
      ReducedImage *ri = *rii;
      TransformedImage *ti = IsTransformedImage(ri);
      if (!ti)
	{
	  GtransfoIdentity identity;
	  ImageGtransfo itf(&identity,&identity, ri->UsablePart(),ri->Name());
	  TIList.push_back(new 
              TransformedImage(TransformedName(ri->Name(),ri->Name()),
			       *ri, &itf)
			   );
	}
      else TIList.push_back(new TransformedImage(*ti));
    }
}

// Name of the transformed with Fakes

static string NameWithFakes(const ReducedImage *Ri)
{
  return ("Fak"+Ri->Name());
}



#include "gtransfo.h"
int Sub::DoOneSub(const ReducedImage *RefStackHere, const ImageSum *NewStackHere, const string &SubName, ImageSubtraction *&SubHere)
{

  SubHere = new ImageSubtraction(SubName, *RefStackHere,*NewStackHere);
  SubHere->FitKernel();
  SubHere->Execute(DoFits | DoCatalog);
  if (AddFakes)
    {
      ReducedImageList tNew = NewStackHere->Components();
      ReducedImageList withoutFakes;
      
      ReducedToTransformed(tNew, withoutFakes);
      ReducedImageList listWithFakes;
      
      double SumNoise=0;
      
      for (ReducedImageIterator tii = withoutFakes.begin(); tii != withoutFakes.end(); ++tii)
	{
	  TransformedImage *ti = dynamic_cast<TransformedImage*>(*tii);
	  Image img = FitsImage((*ti).FitsImageName(Calibrated));
	  Pixel mean, sigma;
	  img.SkyLevel(& mean, &sigma);
	  SumNoise += sigma*sigma;
	  TransformedImageAddFakes *withFakes = new 
	    TransformedImageAddFakes(NameWithFakes(ti),ti->Source(), 
				     ti->Transfo(), FakeList, ti->FromRef()); 
	  listWithFakes.push_back(withFakes);
	}
      // sum up the new stack (with fakes added)
      
      ImageSum newStackWithFakes(NameWithFakes(NewStackHere), listWithFakes);


      

      newStackWithFakes.Execute(DoFits|DoCatalog);
      // borrow the seeing (set when making the catalog)
      //from the new stack w/o fakes:      
      // means that we should put the seeing in the catalog...
      newStackWithFakes.SetSeeing(NewStackHere->Seeing()," picked from the image without fakes");

      FitsImage(newStackWithFakes.FitsImageName(Calibrated)).
	AddOrModKey("TRUESKY",sqrt(SumNoise),
		    " sqrt of the sum of all the images ");
      // The catalog would contain the fakes, not done but anyway useless
      // subtract it using the KernelFit from the regular subtraction.
      ImageSubtraction subtractionWithFakes(NameWithFakes(SubHere), *RefStackHere, newStackWithFakes, SubHere);
      subtractionWithFakes.Execute(DoFits | DoCatalog);
      CandidateStarList foundCands(subtractionWithFakes.AllCandidateCatalogName());
      MatchDetectionsWithFakes( (SEStarList *) &foundCands,
			       subtractionWithFakes.Dir()+"/allcand.matchfakes.list");
      ApplyCutsAndWrite(subtractionWithFakes);
    }

  return 1;
}

// uses SEStarList generically, since we may use this routine both for Candidates and Yquem
void Sub::MatchDetectionsWithFakes(SEStarList *Detections, const string &MatchListName)
{
  double maxdist = 2.0; // pixels
  if (!AddFakes) return;
  SENearStarList putFakes;
  // get the actual fake position/flux in putFakes.
  FakeList->ActualFakes(putFakes); 
  GtransfoIdentity identity;
  StarMatchList *matches = 
          ListMatchCollect(*(BaseStarList *) &putFakes, 
			   *SE2Base(Detections), 
			   &identity, maxdist);

  cout << " for " << MatchListName << " we have " << putFakes.size() 
       << " generated fakes" << endl
       << matches->size() << " detected " << endl;
    
  matches->write(MatchListName);
  delete matches;
}


class CandMatchList : public NStarMatchList {
public:

  // it is extremely important that the lists to match
  //  live at least as long as the NStarMatchList, because the
  // latere deos not hold copies of its inputs. This class serves mainly
  // this purpose.

  vector<CandidateStarList*> candList;
  int nList;

  CandMatchList(const vector<string> &toMatch);

  ~CandMatchList() { for (int i=0; i<nList; ++i) delete candList[i];};


private : 
  // forbid copies, because the default way the compiler does is not gonna work.
  CandMatchList(const CandMatchList &);
  CandMatchList& operator = (const CandMatchList &);

};


CandMatchList::CandMatchList(const vector<string> &toMatch) : NStarMatchList(toMatch.size()), nList(toMatch.size())
{
  DbImage refimage(toMatch[0]);
  
  
  for (int i=0; i<nList ; ++i)
    {
      cout << "Begining the Match of " << toMatch[i] << " with " << toMatch[0] << endl;
      candList.push_back(new CandidateStarList(toMatch[i]));
      BaseStarList* bsl2 = Candidate2Base(candList[i]);
      
            
      // Images are already aligned
      Gtransfo *guess1,*guess2;
      GtransfoIdentity id;
      guess1 = &id;
      guess2 = &id;

      
      SetTransfo(guess1,0,i);
      SetTransfo(guess2,i,0);

      MatchAnotherList(*bsl2,i,2.0,new CandidateStar());
    }
}  

typedef CandMatchList::iterator CandMatchIterator;
typedef CandMatchList::const_iterator CandMatchCIterator;


static double CombinedSTNoise(const NStarMatch &nStarMatch)
{
  double combinedSTNoise = 0;
  for (int i =0; i<nStarMatch.npointers; ++i)
    {
      if ((nStarMatch.StarExist(i)))
	{
	  const CandidateStar *current = (CandidateStar*) nStarMatch.star[i];
	  combinedSTNoise += sqr(current->SigToNoise());
	}
    }
  return sqrt(combinedSTNoise);
}



static  void  MatchThenCutThenWrite(vector<string> toMatch, double STNoiseCut)
{
  CandMatchList cList(toMatch);
  int nList = toMatch.size();
  for (CandMatchIterator si = cList.begin(); si != cList.end(); )
    {
      NStarMatch &sm = *si;
      if (CombinedSTNoise(sm) < STNoiseCut || !sm.StarExist(nList-1))
	si = cList.erase(si);
      else ++si;
    }
  cList.write("AllParameters.list");
  RollingStarList(cList, ReducedImage("ref")).write("Matches.list");

  // Pour Christian : 1 seule detection!
  CandMatchList cListC(toMatch);
  for (CandMatchIterator si = cListC.begin(); si != cListC.end(); )
    {
      NStarMatch &sm = *si;
      if ( CombinedSTNoise(sm) < 5.0 || (sm.actualMatches != 1))
 	si = cList.erase(si);
      else ++si;
    }
  //  cList.write("ChristanAll.list");
  RollingStarList(cListC, ReducedImage("ref")).write("SingleDetections.list");
}


#include "addfakes.h"

int Sub::DoIt()
{
  if (!FixRef)
    FindGeometricReference();  

  ImageSum *RefStack, *NewStack;
  
  //  stack the subparts:
  //int toDo = DoFits | DoCatalog | DoDead | DoSatur;
 
  // The dead is not used anymore (the dead are clipped)
  int toDo = DoFits | DoCatalog | DoSatur | DoDead;

  string NameRef = "ref";
  if (!FileExists(NameRef))
    RefStack = ImagesAlignAndSum(Ref, *GeometricReference, NameRef, toDo);
  else
    RefStack = new ImageSum(NameRef);
    
  // Simulation stuff : computation of the list of fakes
  //  AddFakes=true;
  cout << "AddFakes "<< AddFakes << endl;
  if (AddFakes) 
    {
      FakeList = new SENearStarList;
      
      if (FileExists("infake.list"))
	{
	  cout << "On lit la liste deja existente " << endl;
	  FakeList->read("infake.list");	
	}
      else
	{
	  cout << "On genere La liste de fakes " << endl;
	  if (AssociateGal)
	    MakeListSn(GeometricReference, NewStack, *FakeList, AssociateGal);
	  else
	    MakeListSn(GeometricReference, NewStack, *FakeList);
	  FakeList->write("infake.list");
	}
    }  
      
  // Declaration of the two partial subtractions
  //#ifdef STORAGE
  //ReducedImageList ListOfSub;
  if (!onlyOneSub)
    {
      int n = 1;
      for (ReducedNewIterator it = ListNewList.begin(); it != ListNewList.end();++it )
	{
	  ReducedImageList NewList = *it; 
	  char l[16];
	  sprintf(l, "%d",n);
	  string nameRef = "ref";
	  string nameNew = "new";
	  string nameSub = "sub";
	  nameNew += l;
	  nameSub += l;
	  ImageSum *currentNew;// = NULL;
	  ImageSubtraction *currentSub;
	  if (!FileExists(nameNew))
	    {
	      cout << "Begin the stack of " << nameNew << endl;
	      currentNew = ImagesAlignAndSum(NewList, *GeometricReference, nameNew, toDo);
	    }
	  else
	    {
	      currentNew = new ImageSum(nameNew);
	    }
	  if (!FileExists(nameSub))
	    {
	      cout << "Begin the subtraction " << nameSub << endl;
	      DoOneSub(RefStack,currentNew, nameSub, currentSub);
	      ApplyCutsAndWrite(*currentSub);
	    }
	  else
	    {
	      currentSub = new ImageSubtraction(nameSub,nameRef,nameNew);
	    }
	  ListOfSub.push_back(currentSub);
	  
	  n++;
	  delete currentNew;
	  //delete currentSub;
	}
    }
  // Do the subtraction(s)
  if(!FileExists("new"))
    {
      cout << "Do the stack of all the images " << endl;
      NewStack = ImagesAlignAndSum(AllNew, *GeometricReference, "new", toDo);
    }
  else
    NewStack = new ImageSum("new");
  
  if(!FileExists("sub"))
    {
      DoOneSub(RefStack,NewStack, "sub", Sub);
      ApplyCutsAndWrite(*Sub);
    }
  else 
    {
      string NameRef = "ref";
      string NameNew = "new";
      Sub = new ImageSubtraction("sub",NameRef,NameNew);
    }

  if (NumberOfSub==2)
    ConstructMatchAndCut();

  vector<string> toMatch;
  
  // Do the partial subtraction 
  //#ifdef STORAGE
  if (!onlyOneSub)
    {
      for (ImageSubtractionIterator it = ListOfSub.begin(); it != ListOfSub.end(); ++it)
	{
	  ImageSubtraction *currentSub = *it;
	  string CandCatName = currentSub->AllCandName();
	  if (!FileExists(CandCatName.c_str()))
	    {
	      cerr << "The Candidate Catalog associated to image " << currentSub->Name() << " doesn't exist !! " << endl;
	      if (onlyDet) currentSub->MakeCatalog();
	      else continue;
	    }
	  
	  string DbName = currentSub->Name();
	  toMatch.push_back(CandCatName);
	  // Is deleted there but created two loops before
	  // Due to the constructor of ImageSubtraction
	  delete currentSub;
	}
      //couper
    }
  

  // else ApplyCutsAndWrite(*Sub); // this was I (P.A) think the stuff for 2 subs at the same epoch.
  
/* TO DO 
  - faire la commutation automatique entre mode coincidence ( new1 new2 comme a l'INT )
  et mode rolling
  
  Ca doit marcher en principe, il y a un switch pour le faire si le nombre de soustraction partielle est 2. Julien
  
  - lire le 3.5 qui suit dans les datacards.
  - virer des datacards tout ce qui sert a rien (code + fichier)
  -> J'ai commence a le faire mais pour le faire completement, il faut faire le menage dans les differentes listes.
  
  - virer de YquemStar tout ce qui sert plus a rien.
  -> Voir le point ci-dessus. 
  - est-ce que l'on peut utiliser NStarMatch plutot que YquemStar pour les soustractions traditionnelles.
  
  - et peut etre d'autres choses...

*/


  // this is for the rolling search "
  // 3.5 to be read in the datacards.
  MatchThenCutThenWrite(toMatch,3.5);

  return 1;
}





//pour les cuts
void Sub::ApplyCutsAndWrite(ImageSubtraction &ASubtraction)
{
  YquemStarList stlcand;
  Construction(ASubtraction, stlcand); 
  
  //stlcand.write(ASubtraction.Dir()+"/beforecut.list");
  // Is the same as cand.list
  if (AddFakes)
    MatchDetectionsWithFakes( (SEStarList *) &stlcand,
			      ASubtraction.Dir()+"/cand_cut.matchfakes.list");
  //  string filename = ASubtraction.Dir()+"/candcut.list" ;

  //Apply cuts:
  // on signal to noise as in the datacards (2.5)
  // on the second moment greater than 1e-10
  // on the percentage opf increase (3%)
  // test if is not flagged
  
  DatDetec datdet(DefaultDatacards());
  YquemStarList stlcandcut ;
  stlcand.Cut(stlcandcut,datdet);
  
  // Write
  string CutName = ASubtraction.CandCutName();
  stlcandcut.write(CutName);
  string CutScanName = ASubtraction.CandScanName();
  stlcandcut.write_scan(CutScanName);
  cout << stlcandcut.size() << " candidats selectionnes " << endl ;
  //Cut_Write(stlcand, DirName);
}


// Used for the partial subtraction
void Sub::Cut_Write(CandStarList & stlcand)
{
  DatDetec datdet(DefaultDatacards());
  //datdet.Print();
  
  CandStarList stlcandcut;
  // Cuts : sigtonoise
  //        SigtoNoise1 and SigtoNoise2
  //        saturation on ref
  //        percentage increase
  stlcand.Cut(stlcandcut,datdet);

  cout << stlcandcut.size() << " candidats selectionnes " << endl ;
  stlcandcut.write(Sub->CandCutName());
  stlcandcut.write_scan(Sub->CandCutScanName());
  //  stlcandcut.write_nice(CandCut_NiceList());
  //  stlcandcut.write_scan(CandCut_ScanList());
}



// Associe atous les candidats une etoiles sur la ref
void Sub::Construction(ImageSubtraction &ASub, YquemStarList & stlcand)
{
  DatDetec datdet(DefaultDatacards()); 
  //datdet.Print();

  double seeing_sub = ASub.Seeing();
  double rayon_ref = datdet.n_rayon * seeing_sub; // could be seeing_ref
  
  FitsImage  *pimageref = new FitsImage(ASub.Ref()->FitsName()); // 1 image 
  //cout << ASub.Ref()->FitsName() << endl;
  double backlevel_ref = ASub.Ref()->BackLevel();
  cout << "Pour flux ref : " << rayon_ref << " " << backlevel_ref << endl ;
  
  // Routines in  candstar.cc
  SEStarList refList(ASub.Ref()->CatalogName());
  CandidateStarList stl(ASub.AllCandidateCatalogName());
  
  stlcand.Construct(stl, refList, *pimageref, rayon_ref, backlevel_ref);
  delete pimageref ; 
  stlcand.write(ASub.CandName());
  //  write_scan(stlcand, ASub.CandScanName());
  cout << stlcand.size() << " candidats detectes " << endl ;
}


// Construction + simple
void Sub::ConstructMatchAndCut()
{
  // read the catalog of the sub
  CandidateStarList stl(Sub->AllCandidateCatalogName());

  // reda the catalog of sub1 
  ImageSubtractionIterator it = ListOfSub.begin();
  CandidateStarList stl1((*it)->AllCandidateCatalogName());

  // redad the catalog of sub2
  ++it;  
  CandidateStarList stl2((*it)->AllCandidateCatalogName());  
  

  CandStarList matches;
  
  // Read detection parameters in the datacards
  DatDetec datdet(DefaultDatacards()); 
  //datdet.Print(); Done several times elsewhere

  double seeing_sub = Sub->Seeing();
  double rayon_ref = datdet.n_rayon * seeing_sub; // could be seeing_ref
  
  FitsImage imageref(Sub->Ref()->FitsName()); // 1 image 

  cout << "Pour flux ref : " << rayon_ref << " " << Sub->Ref()->BackLevel()<< endl ;
  SEStarList RefList(Sub->Ref()->CatalogName());
  
  
  matches.Construct_Simple(stl, stl1,  stl2, RefList, imageref, 
			   datdet.dist_min, rayon_ref, Sub->Ref()->BackLevel());

  matches.write(Sub->CandName());

  Cut_Write(matches);


  cout << matches.size() << " candidats detectes " << endl;
}


// ================ DatSim ========================== 

DatSim::DatSim(const string &FileName)
{
  if (FileExists(FileName))
    {
      DataCards cards(FileName);
      LitDataCards(cards);
    }
  else
    {
      cerr << " DatSim::DatSim Cannot open  FileName " << endl;
    }
}

// Reads the parameters for the simulation
void DatSim::LitDataCards(DataCards &Data)
{
  numberOfFakes = Data.IParam("NUMBER_OF_FAKES");
  minMag = Data.DParam("MIN_FAKE_MAG");
  maxMag = Data.DParam("MAX_FAKE_MAG");
  if (minMag > maxMag) swap (minMag,maxMag);
}

void DatSim::Print()
{
  cout << " FAKES : number  of fakes : " <<  numberOfFakes 
      << " mag range : ["<< minMag << ',' << maxMag << ']' << endl;
}


#endif


