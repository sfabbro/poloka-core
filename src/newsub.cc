#include <iostream>
#include <cmath>

#include "imagematch.h"
#include "nstarmatch.h"
#include "newsub.h"
#include "fileutils.h"
#include "quali_box.h"
#include "imagesum.h"
#include "transformedimage.h"
#include "imagesubtraction.h"
#include "fitsimage.h"
#include "allreducedimage.h"
#include "senearstar.h"
#include "transformedimageaddfakes.h"
#include "dimage.h"
#include "toadscards.h"
#include "rollingstar.h"

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
502809o09
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
#502809o09
\endcode

In this specific example, the split of images between NEW1 and NEW2 is
irrelevant since we run a "simple" subtraction (as specified by ONE SUBTRACTION). ADDFAKES is commented and no fakes will be added.


*/


static double sqr(double x){return(x*x);}

// Read the subfile. If you change the syntax here, update the example just above.
NewSub::NewSub(const string &FileName, const bool Overwrite, const bool OnlyDet) : Ref(false), ListNewList(),  AllNew(false), overwrite(Overwrite)
{
  // Boolean for the simulation
  //  FakeList = NULL;

  AssociateGal = false;
  onlyDet = OnlyDet;
  AddFakes = false ;
  FixRef = false;
  FakeList = NULL;   

  //  fit1 = NULL; fit2 = NULL; fitTot = NULL;
  
  FILE *file = fopen(FileName.c_str(),"r");
  if (!file)
    {
      cerr << " Subtraction : cannot open " << FileName << endl;
      return;
    }
  char line[512];
  bool inRef = false;
  bool inNew = false;
  bool inFixRef = false;
  onlyOneSub = false;
  ReducedImageList New(false);
  NumberOfSub = 0;

  while (fgets(line,512,file))
    {

      if (line[0] == '#') continue;
      if (strstr(line,"REF")) 
	{ inFixRef = false; inRef = true; inNew = false; continue;}
      if (strstr(line,"NEW")) 
	{ inFixRef = false; inNew = true; inRef = false; 
	if (NumberOfSub > 0) 
	  { 
	    ListNewList.push_back(New); 
	    New.clear();
	  }
	NumberOfSub++ ;
	continue;
	}
      if (strstr(line,"ONE SUBTRACTION"))
	{ onlyOneSub = true; continue;}
      if (strstr(line,"ADDFAKES")) 
	{ AddFakes = true; cout << "ADDFAKES " << AddFakes << endl; continue;}

      if (strstr(line,"ASSOCIATEGAL"))
	{ AssociateGal = true; continue;}
      if (strstr(line,"FIXGEO"))
	{
	  inNew = false ; inRef = false; inFixRef = true; 
	  FixRef= true; 
	  cout << "The Geometric reference is fixed." << endl; continue;
	}

      if (!inRef && !inNew && !inFixRef)
	{
	  cerr << " ERROR : unrecognised syntax in file " << FileName<< endl ;
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
      if (strlen(start_line) == 0) continue;
  
      ReducedImage *current = new ReducedImage(start_line);
      if (!current->IsValid())
	{
	  cerr << " cannot find DbImage : " << start_line << endl;
          if (onlyDet) { cerr << " tolerable since you asked for detection only" << endl; continue;}
	  exit(1);
	}
      if (inRef)
	{
	  Ref.push_back(current);
	}
      else if (inNew) 
	{
	  New.push_back(current); 
	}
      else if (inFixRef)
	{
	  cout << "Geometric reference will be: " << current->Name() << endl;
	  GeometricReference = current;
	}
      if (inNew)  AllNew.push_back(current->Clone());
      AllImages.push_back(current->Clone());
      
    }
  //cout << AllImages << endl;
  ListNewList.push_back(New); 
  
  if (AllImages.size() == 0)
    {
      cerr << " Did not find any image name in " <<  FileName << " : stop here " << endl;
      exit(1);
    }
  cout << "Beginning the subtraction for "<< NumberOfSub <<" Periods"<< endl;
  fclose(file);
  /* CHECK that all files have a ZEROUSNO fits key */
#ifdef STORAGE
  for (unsigned int i=0 ;i < AllImages.size(); ++i)
    {
      string fitsName = AllImages[i]->FitsName();
      FitsHeader header(fitsName, RW);
      if (AddFakes && !header.HasKey("ZEROUSNO"))
	{
	  //header.AddOrModKey("ZEROUSNO", 33., " bogus zero usno ");
	  cerr << " you requested to add fake SN but there is no ZEROUSNO key in " << fitsName << endl;
	  cerr << " NO FAKES ! " << endl;
	  AddFakes = false;
	}
    }
#endif
  if ( onlyOneSub )
    cerr << "Only One Sub will be processed" << endl ;
}


void NewSub::FindGeometricReference()
{
  
  /* locate the worst seeing image among all */
  GeometricReference = AllImages.front();
  
  for (ReducedImageCIterator ri = AllImages.begin(); ri != AllImages.end(); ++ri)
    {  
      const ReducedImage *current = *ri;
      //cout << current->Name() << endl;
      if (current->Seeing()*current->PixelSize() 
	  > GeometricReference->Seeing()*GeometricReference->PixelSize()) 
	{
	  GeometricReference = current;
	}
    }
  cerr << " Choosing " << GeometricReference->Name() << " as GeometricReference " << endl;
}



static void ReducedToTransformed(const ReducedImageList &RIList, ReducedImageList &TIList)
{
  TIList.clear();
  for (ReducedImageCIterator rii = RIList.begin(); rii != RIList.end(); ++rii)
    {
      const ReducedImage *ri = *rii;
      const TransformedImage *ti = IsTransformedImage(ri);
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
int NewSub::DoOneSub(const ReducedImage *RefStackHere, const ImageSum *NewStackHere, const string &SubName, ImageSubtraction *&SubHere)
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
	  ReducedImage *ri = *tii;
	  TransformedImage *ti = dynamic_cast<TransformedImage*>(ri);
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
void NewSub::MatchDetectionsWithFakes(SEStarList *Detections, const string &MatchListName)
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

int NewSub::DoIt()
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


NewSub::~NewSub()
{
  if (FakeList) delete FakeList;
  int count = 0;
  if (ListOfSub.size() != 0) 
    {
      for (ImageSubtractionIterator it = ListOfSub.begin(); it != ListOfSub.end(); )
	{
	  cout << "~NewSub" << (*it)->Name() << " / count :" << ++count << endl;
	  delete *it;
	  it = ListOfSub.erase(it);	  
	}
    }
}



//pour les cuts
void NewSub::ApplyCutsAndWrite(ImageSubtraction &ASubtraction)
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
void NewSub::Cut_Write(CandStarList & stlcand)
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
void NewSub::Construction(ImageSubtraction &ASub, YquemStarList & stlcand)
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
void NewSub::ConstructMatchAndCut()
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


