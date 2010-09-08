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
#include "swarpstack.h"
#include "imageutils.h"

/*! \page subfile Syntax of the "subfile"
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
is specified by putting "ONE SUBTRACTION" in this file. 
Blank lines in the subfile are ignored.
By default, newmake_sub reads a file named subfile in the current directory.

"FIXGEO" aims at fixing the geometric reference used to align all the images.
Here is an example of a "subfile":

\code
# subtraction for fieldP1ccd9
REF
502809o09 <SubImage>
502816o09
502817o09
NEW1
504604o09
504605o09
NEW2
504613o09
504614o09
ONE SUBTRACTION
#FIXGEO
502809o09
\endcode

In this specific example, the split of images between NEW1 and NEW2 is
irrelevant since we run a "simple" subtraction (as specified by ONE SUBTRACTION). 
Another Example:

 
\code
REF
master_LDP_D3_r_run16 ref=[2262:4544,10686:16051]
NEW SWARP_STACK
744309o00
744309o01
744309o09
744309o10
744309o18
\endcode



*/


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
  while ( (*p != sep))
    {
      if (p == line) return p;
      p--;
    }
  return p+1;
}
 
static string first_word(const char *line, const char *seps = " \t")
{
  int start = strspn(line, seps); // first which is not
  int end = strcspn(line+start,seps);
  return string(line).substr(start,end);;
}
 


// Read the subfile. If you change the syntax here, update the example just above.
Sub::Sub(const string &FileName, const bool Overwrite, const bool OnlyDet) : overwrite(Overwrite)
{
  
  onlyDet = OnlyDet;
  detectOnAllSub = false;
  FixRef = false;
  StringList ToExtract;
  ImageNameToExtract = "";
  onlyOneSub = false;
  
  globnewname = "new";
  globsubname = "sub" ;

  
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
 // check if there is anything left
	  string currentNewName;
	  if (strstr(line,"SWARP_STACK"))
	    {
	      currentNewStack->stackType = SwarpKind;
	      currentNewName = StringToLower(first_word(line));
	    }
	  else currentNewName = StringToLower(last_word(line)); 
	  // remove trailing white space
	  RemovePattern(currentNewName," ");
	  currentNewStack->name = currentNewName;
	  continue;
	}
      if (strstr(line,"ONE SUBTRACTION"))
	{ onlyOneSub = true; continue;}
      if (strstr(line,"DETECT ON ALL SUB"))
	{ detectOnAllSub = true; continue;}
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

      // .. and is not empty

      // check if catalog exists...
      if (current.HasCatalog())
	{
	  SEStarList catalog(current.CatalogName());
	  if (catalog.size() == 0)
	    {
	      std::cerr << " image " << currentName 
			<< " has an empty catalog : drop it" << std::endl;
	      continue;
	    }
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
	{ // search if there is something after the image name
	  std::string remainder = string(start_line + strcspn(start_line," \t"));
	  RemovePattern(remainder," ");
	  if (remainder.size() == 0) Ref.push_back(currentName);
	  else if (strstr(remainder.c_str(), "SubImage") == 0)
	    { // decode name=[a:b,c:d] 
	      char refName[128];
	      int imin,imax,jmin,jmax;
	      /* sscanf format means : everything without "=", 
		 then "=", then you can figure it out */
	      if (sscanf(remainder.c_str(),"%[^=]=[%d:%d,%d:%d]",
			 refName, &imin,&imax,&jmin,&jmax) != 5)
		{
		  std::cerr << " can't decode end of line : " 
			    << remainder << std::endl
			    << line << std::endl
			    << " giving up " << std::endl;
		  exit(1);
		}
	      SubImage ref(refName, currentName, Frame(imin,jmin,imax,jmax));
	      ref.Execute(DoFits | DoWeight | DoSatur | DoCatalog );
	      Ref.push_back(refName);
	      currentName = refName; // for AllInputImages
	      GeomRefName = currentName;
	    }		    
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
      if (AllNew[i].Name() == GlobalNewName())
	{
	  cerr << " \"" << GlobalNewName()<< "\" is a reserved name for the global new image" 
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
	  
static 
void GetAllComponents(vector<NewStack>  & allnew,
			ReducedImageList & allAlignedNew  )
{

  for (unsigned i=0; i< allnew.size(); ++i)
    {
      ReducedImage* stack = allnew[i].newStack; // stack new_i
      ImageSum *sum = dynamic_cast<ImageSum *>(stack);
      if (sum) // this is an actual sum
	{
	  ReducedImageList sumComponents = sum->Components();
	  for (ReducedImageIterator c = sumComponents.begin(); 
	       c != sumComponents.end(); ++c)
	    allAlignedNew.push_back(*c);
	}
      else // the geom ref should be a single image
	{
	      allAlignedNew.push_back(stack);
	}
    }
  return ;
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
	{ string s = *k ; res.push_back(s);}
    }
  cerr << "All New Names : " << res << endl ; 
  return res;
}
	     
// replace image name Original by Substitution in the
// StringList associated to AllNew[i], for all i. and Ref

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
      currentFrame = ApplyTransfo(currentFrame,*current2Large);
      /* if we are at the first image in the loop we have to initialize
	 "toExtract", else we just "increment" it so that we have the union
	 of all input images at the end */
      if (toExtract.Area() == 0.)
	toExtract = (largeFrame*currentFrame);
      else
	toExtract += (largeFrame*currentFrame);
      cout << " after adding " << current.Name() 
	   << ", extract frame = " << toExtract << endl;
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


ReducedImage* Sub::DoOneStack(const StringList &Components, 
			      const string &StackName, 
			      const int ToDo,
			      const StackType ST)
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
  if (ST == RegularKind)
    {
      ReducedImageList unalignedImages(Components);
      ImageSum *stack = ImagesAlignAndSum(unalignedImages, *GeometricReference,
					  StackName, ToDo, RefStack);
      return stack;
    }
  else if (ST == SwarpKind)
    {
      SwarpStack *result = new SwarpStack(StackName, Components, 
					  GeometricReference,
					  GeometricReference->PhysicalSize());
      /* subtractions are usually run in parallel so do not use swarp 
	 internal parallelization :
      */
      AddSwarpDefaultKey("NTHREADS","1");
      result->Execute(ToDo);
      return result;
    }
  else return NULL;
}

static bool FlagDiffractionSpikes(const ReducedImageRef &Ri);


int Sub::DoIt()
{
  FindGeometricReference();
  ReducedImageList NewStacks;
  
  //  stack the subparts:
  //int toDo = DoFits | DoCatalog | DoDead | DoSatur;
 
  // The dead is not used anymore (the dead are clipped)
  int toDo = DoFits | DoCatalog | DoSatur | DoWeight;
  cerr << "Processing ref stack " << endl ;

  if (RefStack ==  NULL) // RefStack isn't alreadt provided.
    {
      RefStack = DoOneStack(Ref, "ref", toDo);
    }
  if (!RefStack)
    {
      cerr << " ERROR: No RefStack, stop here " << endl;
      return 0;
    } 
  //building the stack new_i and the associated sub_i.
  // detection on sub_i
  // aligning in the process the components on ref geom
 
  for (unsigned i=0; i<AllNew.size(); ++i)
    {
      NewStack &stack = AllNew[i]; 
      if (stack.newStack == NULL ) // regular mode, or MC mode in SwarpKind
	{
	  cerr << "Processing " << stack.Name() << endl ;
	  stack.newStack = DoOneStack(stack, stack.Name() , toDo, stack.stackType);
	  FlagDiffractionSpikes(stack.newStack);
	}
      else
	{
	  // MC mode, in RegularKind
	  cerr <<"Processing  in MC mode " << stack.Name() << endl ;
	  if (stack.stackType != RegularKind)
	    cerr <<"Error, bad stackType" << endl ;
	  // pas de SExtractor sur MCnewi
	  (stack.newStack)->Execute(DoFits | DoSatur | DoWeight); 
	  FlagDiffractionSpikes(stack.newStack); //useful?
	}
      
      if ( stack.original_newStack != NULL )
	{
	  // le fits du stack new est fait. on recupere le seeing du stack new d'origine
	  cerr <<"Getting SESEEING value for : " << stack.Name()
	       << " from : " <<  stack.original_newStack->Name() << endl ;
	  string comment = "Seeing from " +  stack.original_newStack->Name() ;
	  stack.newStack->SetSeeing(stack.original_newStack->Seeing(), comment.c_str());
	}
      if ( stack.original_sub == NULL )
	{
	  cerr <<"Building Subtraction " << stack.SubName() << endl ;
	  stack.sub = new ImageSubtraction(stack.SubName(), RefStack, stack.newStack);
	}
      else
	{
	  cerr << "using a previous PsfMatch to build Subtraction " << stack.SubName() <<  endl ;
	  stack.sub = new ImageSubtraction(stack.SubName(), 
					     RefStack, 
					     stack.newStack,
					     stack.original_sub);
	}
      
      stack.sub->Execute(DoFits+DoWeight);
      if (detectOnAllSub) stack.sub->MakeCatalog();
    }
  // building the new stack = stack of all new_i components
  // components already aligned
  // the associated sub and detection.
  if (AllNew.size() > 1)
    {
      ReducedImageList allAlignedNew(false);
      GetAllComponents(AllNew,allAlignedNew);
      GlobalNew = new ImageSum(GlobalNewName(),allAlignedNew, RefStack);     
      if (Original_New != NULL )
	{ 
	  // on recupere le seeing du stack new d'origine
	  GlobalNew->Execute(DoFits); // pas de SExtractor sur MCnew.
	  cerr <<"Getting SESEEING value for : " << GlobalNewName()
	       << " from : " <<  Original_New->Name() << endl ;
	  string comment = "Seeing from " + Original_New->Name() ;
	  GlobalNew->SetSeeing(Original_New->Seeing(), comment.c_str());	  
	}
      else
	 GlobalNew->Execute(DoFits + DoCatalog);

      if (Original_Sub == NULL ) // kernel fit wasn't done before
	GlobalSub = new ImageSubtraction(GlobalSubName(), RefStack, GlobalNew);
      else
	{
	  cerr << "Using a previous PsfMatch to build " 
	       << GlobalSubName() <<  endl ; 
	  GlobalSub = new ImageSubtraction(GlobalSubName(), RefStack, 
					   GlobalNew, Original_Sub);
	}	  
      GlobalSub->Execute(DoFits+DoWeight);
    }
  else // only one subtraction
    {
      // we only have one new/sub
      GlobalNew = AllNew[0].newStack ;
      GlobalSub = AllNew[0].sub;
      if (GlobalSub->Name() != GlobalSubName())
	RenameDbImage(*GlobalSub,GlobalSubName());
      if (AllNew[0].newStack->Name() != GlobalNewName())
	RenameDbImage(*AllNew[0].newStack, GlobalNewName());
    }
  RunDetection();
  return 1;
}

int Sub::ExpectedMagLim()
{
  double nsigma = 3.0;
  size_t Nimages = AllNew.size() + 1;
  static double* photomRatio = new double[Nimages];
  static double* sigmaSeeing = new double[Nimages];
  static double* sigmaSky = new double[Nimages];
  double ZPT;

  ReducedImage ref(Ref.front());
  ReducedImageList unalignedNewImages(&AllNew);

  photomRatio[0] = 1.0;
  sigmaSeeing[0] = ref.Seeing();
  sigmaSky[0] = ref.SigmaBack();
  ZPT = ref.AnyZeroPoint();

  int j=1;
  for (ReducedImageIterator i = unalignedNewImages.begin(); i != unalignedNewImages.end(); ++i)
    {
      ReducedImage *ri = *i;
      double err;      
      photomRatio[j] = QuickPhotomRatio(*ri, ref, err);
      sigmaSeeing[j] = ri->Seeing();
      sigmaSky[j] = ri->SigmaBack();
      j++;
    }

  double sigmaTot2 = 0;
  for (size_t i=0; i<Nimages; ++i)
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
  cerr << "Numbers of Detection on  " <<  GlobalSubName() << " : " << detsOnGlobal.size() << endl ;
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


string NewStack::SubName() const
{ 
  string subname = name ; 
  if (strstr(name.c_str(), "new"))
    return SubstitutePattern(name,"new","sub");
  else
   return(name + "_sub");
}



NewStack::~NewStack()
{
}


/*   Some code related to Diffraction spikes handling 
 It could probably be made clearer and more efficient, but it
moreorless does the job */


class PixCoord 
{
public:
  int x;
  int y;  
  
  PixCoord() { x=y=0;}
  PixCoord(int nx,int ny) {x = nx; y = ny;}

};




class mycluster : public std::list<PixCoord> 
{
public:
  mycluster(int x,int y,Image &image) { addpoint(x,y,image);}
  
  
  void setinimage(Image &image, Pixel value) const 
  {
    	std::list<PixCoord>::const_iterator point = begin();
    	std::list<PixCoord>::const_iterator endpoint = end();
    	for(;point!=endpoint;++point) {image(point->x,point->y)=value;}
  }
  
  
  void enlarge(Image &image) 
  {
    	std::list<PixCoord>::iterator point = begin();
    	std::list<PixCoord>::iterator endpoint = end();
    	for(;point!=endpoint;++point) {addpoint(point->x,point->y,image,false);}   
  }
  
private:
  
  void addpoint(int x,int y,Image &image,bool addthispoint = true) 
  {
    	if(size()>10000)
	{
    	  cerr << "cluster overflow" << endl;
    	  return;
    	}
    	if(image(x,y)==0)return;
    	image(x,y) = 0;
    	if(addthispoint) push_back(PixCoord(x,y));
    	if(y>0) 
	{
      		addpoint(x,y-1,image);
      		if(x>0) addpoint(x-1,y-1,image);
      		if(x<image.Nx()-1)  addpoint(x+1,y-1,image);
    	}
    	if(x>0) addpoint(x-1,y,image);
    	if(x<image.Nx()-1)  addpoint(x+1,y,image);
    	if(y<image.Ny()-1) 
	{
      		addpoint(x,y+1,image);
      		if(x>0) addpoint(x-1,y+1,image);
      		if(x<image.Nx()-1)  addpoint(x+1,y+1,image);
    	}
  }


};

static void findclusters(const Image &image, std::list<mycluster> &clusters) {
  Image newimage = image;
  for (int ix=0;ix< newimage.Nx();ix++)
    for (int iy=0;iy< newimage.Ny();iy++)
      if(newimage(ix,iy)>0) {
	//cout << "New cluster" << endl;
	mycluster cluster(ix,iy,newimage);
	clusters.push_back(cluster);
      }
}


static void saveclustersinimage( const std::list<mycluster> &clusters , Image &image) {
  for (int ix=0;ix< image.Nx();ix++)
     for (int iy=0;iy< image.Ny();iy++)
       image(ix,iy)=0;
  std::list<mycluster>::const_iterator cluster = clusters.begin();
   std::list<mycluster>::const_iterator endcluster = clusters.end();
   for(;cluster!=endcluster;++cluster) {
     cluster->setinimage(image,1);
   }
}

// just a rebinning
static void rebin_image(Image & image, int rebinx, int rebiny) {
  if(rebinx<=0 || rebiny<=0) {
    cerr << "ERROR rebin_image " << rebinx << " " << rebiny << endl;
  }
  // memory consumming 
  int newnx = image.Nx()/rebinx;
  if(image.Nx()%rebinx)
    newnx++;
  int newny = image.Ny()/rebiny;
  if(image.Ny()%rebiny)
    newny++;
  
  Image newimage(newnx,newny);
  Image newweight(newnx,newny);
  Pixel *p = image.begin();
  Pixel *pend = image.end();
  int nx = image.Nx();
  int inx,iny;
  int ix,iy;

  ix=iy=inx=iny=0;
  for (; p!=pend;++p) {
    newimage(inx,iny) += (*p);
    newweight(inx,iny) += 1.;
    ix++;
    if(ix>=nx) {
      ix=0;
      iy++;
      iny = iy/rebiny;
    }
    inx = ix/rebinx;
    //if(inx>= newnx || iny>= newny) {
    // cout << "error " << ix << " " << iy << " = " << inx << " " << iny << endl;
    //}
  }
  newimage/=newweight;
  image=newimage;
} 


// just an unbinning (need something clever here, for the future)
static void unbin_image(Image & image, int nx, int ny, int rebinx, int rebiny) {
  if(nx<=image.Nx() || ny<=image.Ny()) {
    cerr << "ERROR rebin_image " << nx << " " << ny << endl;
  }
  
  // memory consumming 
  Image newimage(nx,ny);
  
  Pixel *p = newimage.begin();
  Pixel *pend = newimage.end();
  int inx,iny;
  int ix,iy;
  
  ix=iy=inx=iny=0;
  for (; p!=pend;++p) {
    (*p) = image(inx,iny);
    ix++;
    if(ix>=nx) {
      ix=0;
      iy++;
      iny = iy/rebiny;
    }
    inx = ix/rebinx;
  }
  image=newimage;
} 





//! enlarge satured clusters in order to get rid of diffraction spikes (aigrettes en francais)
//! method:
//! * rebin satur image (see imagebinning.h)
//! * make clusters 
//! * rebin image (calibrated.fits)
//! * set it to 1 if > N*sigma else 0
//! * enlarge clusters using this new image
//! * save cluster in an image, unbin it back to orginal size and save it as satur.fits.gz 
static bool FlagDiffractionSpikes(const ReducedImageRef &Ri) {
  
  cout << "Entering ReducedImage::FlagDiffractionSpikes" << endl;

  int rebin = 4;
  bool debug = false;
  
  if(!Ri->HasSatur())
    return false;
  
  std::list<mycluster> saturatedclusters;
  
  // read satur and make cluster list
  int nx,ny;
  Image satur; // I want RW access but this is not possible with cfitsio when fitsfile is compressed
  
  {
    FitsImage satur_ro(Ri->FitsSaturName(),RO);
    nx = satur_ro.Nx();
    ny = satur_ro.Ny();
    satur = satur_ro;
  }
  
  if(rebin>1) {      
    rebin_image(satur,rebin,rebin);
  }
   
  // find clusters in saturated image
  findclusters(satur,saturatedclusters);
  //saturatedclusters.sort(&DecreasingClusterSize);
  if(debug) {
    cout << "list of saturated clusters = " << saturatedclusters.size() << endl;
    std::list<mycluster>::iterator cluster = saturatedclusters.begin();
    std::list<mycluster>::iterator endcluster = saturatedclusters.end();
    for(;cluster!=endcluster;++cluster) {
      cout << "size " << cluster->size() << endl;
    }
  }
 
  {
    // copy initial image
    FitsImage initial_image(Ri->FitsName());
    
    Image &newimage =  initial_image;
    
    // now rebin image by a factor rebin
    if(rebin>1) {      
      rebin_image(newimage,rebin,rebin);
    }
    // skylev
    Pixel mean,sigma;
    newimage.SkyLevel(&mean, &sigma);
    // pix = 1 if signal greater than 1.5 sigma
    newimage.Simplify(1.5*sigma,1,0);
    
    // now enlarge satured clusters using this image
    std::list<mycluster>::iterator cluster = saturatedclusters.begin();
    std::list<mycluster>::iterator endcluster = saturatedclusters.end();
    if(debug)
      cout << "list of enlarged clusters" << endl;
    for(;cluster!=endcluster;++cluster) {
      cluster->enlarge(newimage);
      if(debug)
	cout << "size " << cluster->size() << endl;
    } 
  }
  
  // save pixels in clusters in a binary image 
  saveclustersinimage(saturatedclusters,satur);
  
  // un bin image to input size
  if(rebin>1) {   
    unbin_image(satur,nx,ny,rebin,rebin);
  }
  
  // save satur in fits
  {
    FitsHeader shead(Ri->FitsSaturName());
    FitsImage saturfits(Ri->FitsSaturName(),shead,satur); // this is an astuce to write gzip compressed fits image
    saturfits.AddCommentLine("This satur file was modified by FlagDiffractionSpikes");
  }
  cout << "Ending FlagDiffractionSpikes" << endl;
  
  return true;
}
