#include <fstream>
#include "daophot.h"
#include "daophotutils.h"
#include "daophotpsf.h"
#include "daophotio.h"
#include "fileutils.h"

void MakeDaoPsf(ReducedImage &Rim, const bool Redo)
{
  if (FileExists(Rim.ImagePsfName()) && !Redo)
    {
      cout << " MakeDaoPsf : PSF already done for " << Rim.Name() << endl;
      return;
    }

  cout << " MakeDaoPsf : Building PSF" << endl;
  Daophot daophot(Rim);
  daophot.IterPsf(Rim);
  if (!UpdateSeeingFromDaoPsf(Rim)) cerr << " MakeDaoPsf : Image PSF parameters not updated \n";
}


void MakeAllStar(const ReducedImage &Rim, const bool Redo)
{
  string alsfile = Rim.ImageCatalogName(DaophotAls);

  if (FileExists(alsfile) && !Redo)
    {
      cout << " ALLSTAR catalog already done for " << Rim.Name() << endl;
      return;
    }

  if (!FileExists(Rim.ImagePsfName()))
    {
      cout << " ALLSTAR : DAOPHOT PSF not yet done for " << Rim.Name() << endl;
      return;
    }

  cout << " Producing ALLSTAR catalog" << endl;
  Daophot daophot(Rim);
  daophot.AllStar(Rim.CatalogName());
  cout << " ALLSTAR complete " << endl;
}


void MakeDaoPsfCat(ReducedImage &Rim, const bool DoPsf, const bool DoCat,  
		   const bool Manually, const bool Combine, const bool Redo)
{
  bool dopsf = DoPsf;
  if (DoPsf && FileExists(Rim.ImagePsfName()) && !Redo)
    {
      cout << " MakeDaoPsfCat : PSF already done for " << Rim.Name() << " : only updating." << endl;
      if (!UpdateSeeingFromDaoPsf(Rim)) 
	cerr << " MakeDaoPsfCat : Image PSF parameters not updated \n";
      dopsf = false;
    }
  
  bool docat = DoCat;
  if (DoCat && FileExists(Rim.ImageCatalogName(DaophotAls)) && 
      FileExists(Rim.ImagePsfName()) &&!Redo)
    {
      cout << " ALLSTAR catalog already done for " << Rim.Name() << endl;
      docat = false;
    }

  if ((dopsf) || (docat))
    {
      SEStarList sex(Rim.CatalogName());
      write_dao<DaophotAp>(Rim, sex);
      Daophot daophot(Rim);
      if (dopsf) 
	{
	  daophot.IterPsf(Rim, Manually);
	  if (!UpdateSeeingFromDaoPsf(Rim)) 
	    cerr << " MakeDaoPsfCat : Image PSF parameters not updated \n";
	}      
      if (docat) daophot.AllStar(Rim.ImageCatalogName(DaophotAp));
      cout << " DAOPHOT PSF and ALLSTAR catalog complete " << endl;
    }

  if (Combine && FileExists(Rim.ImageCatalogName(DaophotAls))) 
      CombineSextractorDaophot<DaophotAls>(Rim);
}

void MakePrecisePsf(const ReducedImage &Rim)
{
  cout << " MakePrecisePsf : Building PSF for " << Rim.Name() << endl;
  Daophot daophot(Rim);
  daophot.PrecisePsf(Rim);
}

