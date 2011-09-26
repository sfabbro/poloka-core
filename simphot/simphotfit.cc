#include "simphotfit.h"
#include "imagepsf.h"
#include "lightcurvesyntaxerror.h"

  
#include <fstream>
#include <iomanip>
#include <string>

/* extracted from SimPhotFit constructor, so that it can be called
   when parsing the lightcurve file. It enables to check values when
   parsing the input file.   
*/

int DecodeFitType(const int FitType)
{
  int toDo = FIT_SKY;
  switch (FitType)
    {
    case -2 : toDo = FIT_FLUX + FIT_POS; break; // NO SKY
    case -1 : toDo += (FIT_FLUX + FIT_GALAXY); break; 
    case 0 : toDo += (FIT_FLUX + FIT_GALAXY + FIT_POS); break; 
    case 1 : toDo += (FIT_FLUX + FIT_POS); break; 
    case 2 : toDo += (FIT_GALAXY); break;
    case 3 : toDo += FIT_FLUX; break;
    default : 
      {
	char mess[128];
	sprintf(mess,"DecodeFitType : Wrong FitType provided (%d) ",FitType);
	cout << mess << endl;
	cout << " available FitTypes : "<< endl;
	cout << " -1 : fit galaxy, fluxes at fixed (provided) position " << endl;
	cout << " 0 : fit galaxy, fluxes, and position" << endl;
	cout << " 1 : fit fluxes and position, without galaxy " << endl;
	cout << " 2 : fit galaxy only " << endl;
	cout << " 3 : fit fluxes, at fixed position, w/o galaxy " << endl;
	throw (LightCurveSyntaxError("",mess));
      }
    }
  return toDo;
}


static string todo_2_string(const int ToDo)
{
  char s[64];
  sprintf(s,"g%d s%d f%d p%d",
	  (ToDo & FIT_GALAXY)>0, 
	  (ToDo &FIT_SKY)>0,
	  (ToDo & FIT_FLUX)>0, 
	  (ToDo & FIT_POS)>0);
  return string(s);
}

/*************************** SimPhotFit class routines *****************/

SimPhotFit::SimPhotFit(const ObjectToFit &Obj, const LightCurveFile &LCFile):
  Model(LCFile, Obj),
  objToFit(Obj), lcFile(LCFile)
{ 
  vignette_size_n_seeing = 4. ;// vignette size in seeing

  maxKernelSize = 0;
  toDo = DecodeFitType(objToFit.FitType());
}


bool SimPhotFit::DoTheFit()
{
  // determine the size of the model we are going to fit
  // we have to load the PSF of the geomRef at the location of the SN
  BuildVignettes();
  //  cout << "BuildVignettes OK " << endl;
  if (toDo & FIT_GALAXY) FindModelBoundaries();
  //  cout << "FindModelBoundaries OK " << endl;
  UpdateResiduals(); // should become useless
  //  cout << "UpdateResiduals OK " << endl;
  bool ok = true;
  // one can only fit position if fluxes are somehow defined:
  // so we first don't fit position to get fluxes

  if (toDo & FIT_POS)
    {
      cout << "DoTheFit : First minimization at fixed position" << endl;
      ok = OneMinimization(toDo & ~FIT_POS & FIT_ALL,1,1.); 
    }

  cout << "DoTheFit : do the requested fit : " << todo_2_string(toDo) << endl;
  ok = (ok && OneMinimization(toDo,15,0.1));
  cout << "OneMinization: " << toDo << " " << ok << endl;
  if (!ok) return false;

  //TODO : put 5 in the datacards
  double nsig = 5;
  cout << "DoTheFit :  removing outliers @ " << nsig << " sigmas " << endl;
  int outPix = 0;
  for (VignetteIterator i = vignetteList.begin(); i != vignetteList.end(); ++i)
    outPix += (*i)->KillOutliers(nsig);
  cout << " DoTheFit : removed " << outPix << " pixels in total " << endl; 

  // we have to refit, either if we discarded pixels or if the fit 
  // is non linear (earlier convergence is crude)
  if ( (toDo & FIT_POS) || outPix != 0) 
    ok = (ok && OneMinimization(toDo,15,0.01));

  // A has to be inverted at the end of the fit to get variances
  cout << " inverting matrix " << endl;
  cholesky_invert(A,"U");
  
  return ok;
}


bool SimPhotFit::BuildVignettes()
{
  const RImageList &images = lcFile.Images();
  bool willFitGal = (toDo & FIT_GALAXY);
  for (RImageCIterator i = images.begin(); i != images.end(); ++i)
    {
      const RImageRef &current = *i;
      Vignette *v = new Vignette(*this, current);
      // we need the convolution kernel only if we fit the galaxy
      if (v->SetGeomTransfos() && (!willFitGal || v->SetKernel()))
	{
	  vignetteList.push_back(v);
	  maxKernelSize = max(maxKernelSize,v->HalfKernelSize());
	}
      else
	{
	  cout << " dropping image " << current->Name() << endl;
	  delete v;
	}
    }
  vignetteList.sort(IncreasingDate);
  return true;
}

#include "matvect.h"

bool SimPhotFit::OneMinimization(const int CurrentToDo, const int MaxIter, 
				 const double DeltaChi2)
{
  double oldChi2; int oldDof;
  CumulateChi2(oldChi2,oldDof);
  double chi2; int ndof;
  cout << "OneMinimization CumulateChi2 ok, oldChi2, oldDof " << setprecision(15) << oldChi2 << " " << oldDof << endl;
  int niter = 0;
  while (niter<MaxIter)
    {
      AssignIndicesAndToDo(CurrentToDo);
      if (posIndex != -1 && fluxMap.size() == 0)
	{
	  cout << "SimPhotFit::OneMinimization : you are trying to fit " << endl
	       << "an object position and all epochs are assumed to have a null flux ..." << endl
	       << " The fit aborts here. It could be that MJD are messed up " << endl;
	  cout << " try header -k TOADMJD <fitsfile> on a sample of input images" << endl;
	  return false;
	}
      if (!OneIteration(CurrentToDo))  return false;
      CumulateChi2(chi2, ndof, true );
      niter++;

      cout << "OneMinimization : niter oldchi2 chi2 " << niter << setprecision(15) << " " << oldChi2 << " " << chi2 << endl ;
      if (chi2 < oldChi2-DeltaChi2)
	{
	  cout << "OneMinimization : Chi2 decreased (chi2<oldChi2-DeltaChi2) : delta chi2 =" << oldChi2-chi2 << endl;
	  oldChi2 = chi2; 
	  continue ; 
	}
      else // chi2 increased or is slightly decreased
	{
	  cout << "OneMinimization : chi2>oldChi2-DeltaChi2 --> parabolic approx : delta chi2 =" << oldChi2-chi2 << endl;
	  double c0=oldChi2,c1,c2=chi2;
	  double x0=0,x1=0.5,x2=1.;
	  DispatchOffsets(B,-0.5, false);
	  CumulateChi2(c1,ndof,false);

	  for (int i=0; i < 10; ++i)
	    {
	      cout << "OneMinimization : c0, c1, c2 " << c0 << " " << c1 << " " << c2 << endl ;
	      //Need a break criteria if approximation had converged
	      if ( abs(c0-c2)<DeltaChi2 )
		{
		  if (abs(c0-c1)<DeltaChi2)
		    {
		      cout << "OneMinimization : Parabolic approx. converged (|c0-c2|<DeltaChi2 && |c0-c1|<DeltaChi2) :  c0-c2=" 
			   << c0-c2 << " c0-c1=" << c0-c1 << endl;
		      break;
		    }
		  else // au cas ou c0 et c2 sont de chaque cote du puits de parabole, a la meme hauteur
		    cout << "OneMinimization  : Parabolic approx. near convergence (|c0-c2|<DeltaChi2) :  c0-c2=" 
			   << c0-c2 << " c0-c1=" << c0-c1 << endl;
		}
	      //Assume Chi2 as a parabole y = a*x^2 + b*x + c
	      //Then xmin = -b/(2*a)
	      double a = (c0*(x1-x2) - c1*(x0-x2) + c2*(x0-x1)) / ((x0*x0 - x0*(x1+x2) + x1*x2) * (x1-x2));
	      double b =  -(c0*(x1*x1-x2*x2) - c1*(x0*x0-x2*x2) + c2*(x0*x0-x1*x1)) / ((x0*x0 - x0*(x1+x2) + x1*x2) * (x1-x2));
	      double c = (c0*x1*x2*(x1-x2) - x0*(c1*x2*(x0-x1) - c2*x0*(x0-x1)))   / ((x0*x0 - x0*(x1+x2) + x1*x2) * (x1-x2));
	      if (fabs(a) < 1e-20)
		{
		  cout << "OneMinimization  : zero value for a : should be protected ????? " << endl ;
		}
	      
	      if ( a < 0 )
		{
		cout << "OneMinimization : negative value for a : c0= " 
		     << c0 << " c1= " << c1  << " c2= " << c2 << endl ;
		}

	      double xmin =  -b/(2*a);
	      printf("OneMinization : Coef. of parabolic approximation of Chi2 : a=%E b=%E c=%E \n",a,b,c);
	      printf("OneMinization : x0 x1 x2 xmin %E %E %E %E \n", x0,x1,x2,xmin);
	      printf("OneMinization : c0 c1 c2 %E %E %E\n",c0,c1,c2);

	      if (xmin<x0)
		{
		  x1 = x0;
		  c1 = c0;
		  DispatchOffsets(B,xmin-x1, false);
		  CumulateChi2(c0,ndof,false);
		  x0 = xmin;
		}
	      else if (xmin>x2)
		{
		  x1 = x2;
		  c1 = c2;
		  DispatchOffsets(B,xmin-x1);
		  CumulateChi2(c2,ndof,false);
		  x2 = xmin;
		}
	      else 
		{
		  if (xmin<x1)
		    {
		      x2 = x1; c2 = c1;
		    }
		  else
		    {
		      x0 = x1; c0 = c1;
		    }
		  DispatchOffsets(B,xmin-x1, false);
		  CumulateChi2(c1,ndof,false);
		  x1 = xmin;
		}
	      chi2 = c1;
	      cout << "OneMinization : parabolic approx :  chi2,  step "<< chi2 << ' ' << x1 << endl;

	    }//end loop on straight line optimization

	  if (chi2<=oldChi2)
	    {
	      if (oldChi2-chi2 < DeltaChi2) 
		{
		  cout << " OneMinization : converged : niter oldchi2, chi2, delta chi2 =" << niter << " " << oldChi2<< " " << chi2 << " " << oldChi2-chi2 << endl; 
		  break;
		}
	    }

	  oldChi2 = chi2;
	
	}// end else chi2 has slightly decreased or increased

    }//end loop while

  return true;
}

bool SimPhotFit::OneIteration(const int CurrentToDo)
{
  A.allocate(matSize,matSize);
  B.allocate(matSize);
  for (VignetteIterator i = vignetteList.begin(); i != vignetteList.end();++i)
    {
      Vignette &v = **i;
      if (toDoMap.find(&v) == toDoMap.end())
	{
	  cout << " BIZARRE : vignette dans la liste et pas dans la map..." << endl;
	  continue;
	}
      int toDo = toDoMap[&v];
      //      cout << v.Name() << ": " << todo_2_string(toDo) << endl;
      v.FillAAndB(A,B,toDoMap[&v]);
    }

  if (!Solve(A,B,"U", (CurrentToDo & FIT_GALAXY), maxKernelSize)) return false;

  // push results into the right places
  DispatchOffsets(B);
  return true;
}

void SimPhotFit::UpdateResiduals()
{
  for (VignetteIterator i = vignetteList.begin(); i != vignetteList.end();++i)
    {
      Vignette &v = **i;      
      v.UpdateResiduals();
    }
}

void SimPhotFit::CumulateChi2(double &Chi2, int &NDof, const bool Print) const
{
  Chi2 = 0;
  NDof = -matSize;
  for (VignetteCIterator i = vignetteList.begin(); i != vignetteList.end();++i)
    {
      const Vignette &v = **i;
      Chi2 += v.Chi2();
      NDof += v.NTerms();
    }
  if (Print)
    cout << " **** chi2/ndof = " << Chi2 << '/' << NDof << " = " << Chi2/NDof << endl;
}





static int nearest_integer(const double &x)
{
  return int(floor(x+0.5));
}

#include "frame.h"

void SimPhotFit::FindModelBoundaries()
{
  Frame modelFrame;
  for (VignetteCIterator i = vignetteList.begin(); i != vignetteList.end();++i)
    {
      const Vignette &v = **i;
      Frame frame;
      v.ComputeModelLimits(frame);
      //      cout << "# FindModelBoundaries - vignette frame size : " << frame.Nx() << endl ;
      modelFrame += frame;
    }
  cout << " Limits (in pixels) of the fitted model " << modelFrame;
  galaxyPixels.Allocate(nearest_integer(modelFrame.xMin), 
			nearest_integer(modelFrame.yMin), 
			nearest_integer(modelFrame.xMax),
			nearest_integer(modelFrame.yMax));
  galaxyPixels.SetVal(0.);
}




// what concerns indices in A and B 

/* assignment routine : ordering : galaxy, sky ('s), flux ('),
position.  If one set is not fitted, it just disappears.  Do not
change the ordering: it is used to fill only one half of the
matrix. (Vignette::FillAandB).
*/
void SimPhotFit::AssignIndicesAndToDo(const int CurrentToDo)
{
  skyMap.clear();
  fluxMap.clear();
  toDoMap.clear();
  posIndex = -1;
  nParamGal = 0;
  int freeIndex = 0;
  // galaxy
  if (CurrentToDo & FIT_GALAXY) 
    {
      nParamGal = galaxyPixels.Ntot();
      freeIndex += nParamGal;
    }
  // sky 's
  // one sky has to be fixed, or the problem is degenerate.
  bool fixedOneSky = false;
  for (VignetteIterator i = vignetteList.begin(); i != vignetteList.end();++i)
    {
      Vignette &v = **i;
      int toDoVignette = CurrentToDo & v.CanDo();

      //cout << "toDoVignette CurrentToDo V.CanDo" << toDoVignette << CurrentToDo << v.CanDo() << endl ;
  
      if (toDoVignette & FIT_SKY)
	{
	  if (!fixedOneSky)
	    {
	      toDoVignette &= (~FIT_SKY & FIT_ALL);
	      fixedOneSky = true;
	    }
	  else
	    skyMap[&v] = freeIndex++;
	}
      toDoMap[&v] = toDoVignette;
      //cout << "toDoMap[&v]" << toDoMap[&v] << endl ;
    }  

  // fluxes
  for (VignetteIterator i = vignetteList.begin(); i != vignetteList.end();++i)
    {
      Vignette &v = **i;
      if (toDoMap.find(&v) == toDoMap.end()) abort(); // should never happen
      int toDoVignette = toDoMap[&v];
      //cout << "toDoVignette :" << toDoVignette << endl ;
      if (toDoVignette & FIT_FLUX) fluxMap[&v] = freeIndex++;
      //cout << "toDoVignette FIT_FLUX fluxMap[&v] :" << toDoVignette << FIT_FLUX << fluxMap[&v] << endl ;
    }
  // position
  if (CurrentToDo & FIT_POS) 
    {
      posIndex = freeIndex;
      freeIndex +=2;
    }
  matSize = freeIndex;
  cout << " SimPhotFit::AssignIndicedAndToDo: request = " << todo_2_string(CurrentToDo) << endl;
  cout << " SimPhotFit::AssignIndicesAndToDo: number of parameters : " 
       << matSize
       << " g:" << nParamGal 
       << " s:" << skyMap.size()
       << " f:" << fluxMap.size()
       << " p:" << int((posIndex == -1)? 0 : 2) 
       << endl;

}

int SimPhotFit::SkyIndex(const Vignette* V) const
{
  VignetteMap::const_iterator i = skyMap.find(V);
  return (i==skyMap.end()) ?  -1 : i->second;
}

int SimPhotFit::FluxIndex(const Vignette* V) const
{
  VignetteMap::const_iterator i = fluxMap.find(V);
  return (i==fluxMap.end()) ?  -1 : i->second;
}


void SimPhotFit::DispatchOffsets(const Vect& Offsets, const double Fact,
				 const bool Verbose)
{
  if (nParamGal)
    {
      if (nParamGal != galaxyPixels.Ntot())
	{
	  cout << " inconsistency in SimPhotFit::DispatchOffsets" << endl;
	  abort();
	}
      int k = 0;
      for (int j=galaxyPixels.ymin; j <galaxyPixels.ymax; ++j)
	for (int i=galaxyPixels.xmin; i <galaxyPixels.xmax; ++i)
	  galaxyPixels(i,j) += Fact*Offsets(k++);
      // we now have an actual galaxy
      hasGalaxy = true;
    }

  // update skies
  for (VignetteMap::const_iterator i = skyMap.begin(); i != skyMap.end(); ++i)
    {
      // cast violates constness...
      Vignette *v = (Vignette *)i->first;
      v->SetSky(v->GetSky()+ Fact*Offsets(i->second));
    }

  //update fluxes. 
  // Use the VignetteList loop to get the printouts in the list order
    for (VignetteIterator i = vignetteList.begin(); i != vignetteList.end();++i)
      {
      Vignette *v = *i;
      int fluxIndex = FluxIndex(v);
      if (fluxIndex < 0) continue;
      v->SetFlux(v->GetFlux()+ Fact * Offsets(fluxIndex));
      if (Verbose)
	cout << " current : " << v->Name() << " f= " << v->GetFlux() 
	     << " s = " << v->GetSky() << endl;
    }
  if (posIndex >= 0)
    {
      Point oldPos = objectPos;
      Point delta(Fact * Offsets(posIndex), Fact * Offsets(posIndex+1));
      objectPos = oldPos + delta;
      if (Verbose)
	cout << " new position : " << ObjectPosInImage() << " delta " << delta << endl;
    }
  UpdateResiduals();
}


/*** IO's of simultaneous fit *******
     We have a "detailed" IO mode for SNe and other single objects
     (SimPhotFit::Write()),
     and the calibration IO mode, for which we build an ntuple
     SimPhotFit::WriteTupleHeader and SimPhotFit::WriteTupleEntries()).
     The calls to one or the other is in LightCurveFile class

*/



#define FLUX_WEIGHT_NAME "lightcurve.weight.dat"


bool SimPhotFit::Write(const string &Dir, const bool WriteVignettes, const bool WriteMatrices)
{

  unsigned nflux = fluxMap.size();
  VignetteMap newFluxMap;
  unsigned *indexMapping = new unsigned[nflux];
  Point posInImage=ObjectPosInImage();

  ofstream lc((Dir+"lc2fit.dat").c_str());
  lc << "@INSTRUMENT MEGACAM" << endl;
  lc << "@BAND " << RefImage().Band() << endl;
  lc << "@MAGSYS AB" << endl;
  lc << "@COORD_X " << posInImage.x << endl;
  lc << "@COORD_Y " << posInImage.y << endl;
  lc << "@REFERENCEIMAGE " << RefImage()->Name() << endl;
  lc << "@WEIGHTMAT " << FLUX_WEIGHT_NAME << endl;
  lc << "@PSFZPERROR 12." << endl;
  lc << "# Date :"  << endl;
  lc << "# Flux :"  << endl;
  lc << "# Fluxerr :"  << endl;
  lc << "# ZP:"  << endl;
  lc << "# sky :" << endl;
  lc << "# esky :" << endl;
  lc << "# seeing: " << endl;
  lc << "# exptime: " << endl;
  lc << "# phratio: photom ratio" << endl;
  lc << "# gseeing: GFSeeing" << endl;
  lc << "# sesky: SESky" << endl;
  lc << "# sigsky: SIGSky" << endl;
  lc << "# name :" << endl;
  lc << "#end" << endl;

  //  make up a zp
  double zp = 12;
  // don't use this loop here because the ordering is not controlled
  //  for (VignetteMap::const_iterator i = fluxMap.begin(); 
  //   i != fluxMap.end(); ++i, ++k)
  unsigned k = 0;

  if ( vignetteList.size() != nflux ) cerr << "IndexMapping(vignetteList.size())  nflux(fluxMap.size()) " <<  vignetteList.size() << " " << nflux << endl ;

  for (VignetteCIterator i = vignetteList.begin(); i != vignetteList.end(); ++i)
    {
      const Vignette *v = *i;

      string dbim_name= v->Name();

      int fluxIndex = FluxIndex(v);
      if (fluxIndex <0 )  continue;
      double sigSky = 0;
      int skyIndex = SkyIndex(v);
      if (skyIndex>=0) sigSky = sqrt(A(skyIndex,skyIndex));
      lc << setw(14) << setprecision(11) << v->MJD() << ' ' 
	 << setw(14) << setprecision(9) << v->GetFlux() << ' ' 
	 << setw(12) << setprecision(9) << sqrt(A(fluxIndex,fluxIndex)) << ' '
	 << setw(7) << setprecision(5) << zp << ' '
	 << setw(10) << setprecision(7) << v->GetSky() << ' '
	 << setw(10) << setprecision(7) << sigSky << ' '
	 << setw(8) << setprecision(7) << v->Seeing() << ' '
	 << setw(8) << setprecision(7) << v->ExpTime() << ' '
	 << setw(8) << setprecision(7) << v->PhotomRatio() << ' '
	 << setw(8) << setprecision(7) << v->GFSeeing() << ' '
	 << setw(8) << setprecision(7) << v->SESky() << ' '
	 << setw(8) << setprecision(7) << v->SIGSky() << ' '
	 << dbim_name << ' '
	 << endl;
      newFluxMap[v] = k;
      indexMapping[k] = fluxIndex;
      k++;
  }
  lc.close();

  ofstream lstream((Dir+"lightcurve_sn.dat").c_str());
  lstream << "# x : x position (pixels) "  << endl;
  lstream << "# y : y position (pixels) "  << endl;
  lstream << "# flux : flux en unites du pixel "  << endl;
  lstream << "# sky : mean sky value per pixel "  << endl;
  lstream << "# varflux : variance of measured flux " << endl;
  lstream << "# varx : variance in x position (pixels) " << endl;
  lstream << "# vary : variance in y position (pixels) " << endl;
  lstream << "# varsky : variance in sky measurement " << endl;
  lstream << "# mjd : modified julian date of the observation" << endl; 
  lstream << "#end" << endl;

  for (VignetteCIterator i = vignetteList.begin(); i != vignetteList.end(); ++i)
    {
      
      const Vignette *v = *i;

      double sigPosX = 0, sigPosY = 0;
      int posIndex = PosIndex();
      if (posIndex >= 0)
	{
	  sigPosX = sqrt(A(posIndex,posIndex));
	  sigPosY = sqrt(A(posIndex+1,posIndex+1));
	}

      int fluxIndex = FluxIndex(v);
      if (fluxIndex < 0) continue;
      double sigFlux=sqrt(A(fluxIndex,fluxIndex));

      double sigSky = 0;
      int skyIndex = SkyIndex(v);
      if (skyIndex >= 0) sigSky = sqrt(A(skyIndex,skyIndex));

      Point fittedPos( ObjectPosInImage());

      lstream << setw(10) << setprecision(7) << fittedPos.x << ' ' 
	      << setw(10) << setprecision(7) << fittedPos.y << ' '
	      << setw(10) << setprecision(7) << v->GetFlux() << ' ' 
	      << setw(10) << setprecision(7) << v->GetSky() << ' ' 
	      << setw(10) << setprecision(7) << sigFlux << ' '
	      << setw(10) << setprecision(7) << sigPosX  << ' ' 
	      << setw(10) << setprecision(7) << sigPosY << ' '
	      << setw(10) << setprecision(7) << sigSky << ' '
	      << setw(13) << setprecision(10) << v->MJD() << ' '

	      << endl;
    }

  lstream.close();

  //flux covariance matrix
  Mat fluxWeight(nflux,nflux);
  // A is a covariance, so fluxWeight is also a covariance 
  A.ExtractSubMat(indexMapping, nflux, fluxWeight);
  fluxWeight.writeFits("flux_pmat_sn.fits");


  // convert it to weight
  fluxWeight.CholeskyInvert("U");
  fluxWeight.Symmetrize("U");
  
  fluxWeight.writeASCII(Dir+FLUX_WEIGHT_NAME);

  delete [] indexMapping;

  ofstream vi((Dir+"vignetteinfo.dat").c_str());
  vi << "# Date : (days since Jan  1st, 2003)"  << endl;
  vi << "# Flux :"  << endl;
  vi << "# Fluxerr :"  << endl;
  vi << "# sky :" << endl;
  vi << "# esky : " << endl;
  vi << "# seeing : pixels sigma" << endl;
  vi << "# chi2 : contribution to the total chi2" << endl;
  vi << "# ndof : number of pixels of the stamp" << endl;
  vi << "# exptime : in seconds" << endl;
  vi << "# fitgal : " << endl;
  vi << "# fitsky : " << endl;
  vi << "# fitflux : " << endl;
  vi << "# fitpos : " << endl;
  vi << "# indexlc : position in the lc and covariance matrix files" << endl;
  vi << "# name : image name" << endl;
  vi << "#end" << endl;

  for (VignetteCIterator i = vignetteList.begin(); i != vignetteList.end();++i)
    {
      const Vignette *v = *i;

      string dbim_name= v->Name();

      int flag = 0; 
      if (toDoMap.find(v) != toDoMap.end()) flag = toDoMap[v];

      double sigFlux = 0;
      int fluxIndex = FluxIndex(v);
      if (fluxIndex >= 0)  sigFlux=sqrt(A(fluxIndex,fluxIndex));

      double sigSky = 0;
      int skyIndex = SkyIndex(v);
      if (skyIndex >= 0) sigSky = sqrt(A(skyIndex,skyIndex));

      int indexLC = -1;
      if (newFluxMap.find(v) != newFluxMap.end()) 
	indexLC = newFluxMap[v];
      vi << setw(10) << setprecision(8) << v->MJD() << ' ' 
	 << setw(11) << setprecision(9)<< v->GetFlux() << ' ' 
	 << setw(11) << setprecision(9) <<sigFlux << ' ' 
	 << setw(10) << setprecision(6) <<v->GetSky() << ' '
	 << setw(10) << setprecision(6) <<sigSky << ' '
	 << setw(6) << setprecision(4) <<v->Seeing() << ' '
	 << setw(8) << setprecision(6) <<v->Chi2() << ' '
	 << setw(8) << setprecision(5) << v->NTerms() << ' ' 
	 << setw(8) << setprecision(6) <<v->ExpTime() << ' '
	 << ((flag&FIT_GALAXY)>0) << ' '
	 << ((flag&FIT_SKY)>0) << ' '
	 << ((flag&FIT_FLUX)>0) << ' '
	 << ((flag&FIT_POS)>0) << ' '
	 << indexLC << ' '
	 << dbim_name << ' '
	 << endl;
    }
  vi.close();



  //galaxy vignette and eventually residual vignettes
  if (toDo & FIT_GALAXY)
    {
      galaxyPixels.WriteFits(Dir+"galaxy_sn.fits");
      PixelBlock galaxyWeight(galaxyPixels);
      PIXEL_LOOP(galaxyWeight,i,j)
      {
	int pixelIndex = galaxyPixels.PixelIndex(i,j);
	galaxyWeight(i,j) = 1./A(pixelIndex,pixelIndex);
      }
      galaxyWeight.WriteFits(Dir+"galaxy.weight.fits");
    }
  if (WriteVignettes)
    for (VignetteCIterator i = vignetteList.begin(); i != vignetteList.end();++i)
      {	  
	const Vignette *v = *i;
	if (fluxMap.find(v) != fluxMap.end())
	  {
	   const Vignette &v2 = **i;
	    v2.Write(Dir);
	  }
      }

  if (WriteMatrices)
    {

      // get covariance matrix
      A.writeFits(Dir+"pmat_sn_bkup.fits");

      // create and get vector of flux
      int i=0;
      for (VignetteCIterator it = vignetteList.begin(); it != vignetteList.end() ; ++it)	
	{   
	  const Vignette *v = *it;
	  if (fluxMap.find(v) != fluxMap.end())
	    {
	      if ((toDo&FIT_FLUX)&&(v->mightFitFlux))
		i++;
	    }
	}

      Mat vm(1,i);

      i=0;
      for (VignetteCIterator it = vignetteList.begin(); it != vignetteList.end() ; ++it)
	{
	  const Vignette *v = *it;
	  if (fluxMap.find(v) != fluxMap.end())
	    {
	      if ((toDo&FIT_FLUX)&&(v->mightFitFlux))
		{
		  vm(0,i)= v->GetFlux();
		  cout << "Flux=" << v->GetFlux() << endl; 
		  i++;
		}
	    }
	}

      vm.writeFits(Dir+"/vec_sn.fits");


      // create and get matrix of Nights and Images
      FillNightMat().writeFits(Dir+"/nightmat_sn.fits");
    }

  return true;
}


void SimPhotFit::WriteTupleHeader(ostream &Stream, const int NStars) const
{
  Stream << "@NSTARS " << NStars << endl;
  Stream << "@NIMAGES " << vignetteList.size() << endl;
  Stream << "#x :" << endl;
  Stream << "#y :" << endl;
  Stream << "#xerror :" << endl;
  Stream << "#yerror :" << endl;
  Stream << "#flux :" << endl;
  Stream << "#error :" << endl;
  Stream << "#sky :" << endl;
  Stream << "#skyerror :" << endl;
  Stream << "#mjd : obs date " << endl;
  Stream << "#seeing: " << endl;
  Stream << "#exptime: " << endl;
  Stream << "#phratio: photom ratio" << endl;
  Stream << "#gseeing: GFSeeing" << endl;
  Stream << "#sesky: SESky" << endl;
  Stream << "#sigsky: SIGSky" << endl;
  Stream << "#mag :" << endl;
  Stream << "#mage :" << endl;
  Stream << "#ra : initial " << endl;
  Stream << "#dec : initial " << endl;
  Stream << "#ix : initial x" << endl;
  Stream << "#iy : initial y" << endl;
  Stream << "#u : from catalog" << endl;
  Stream << "#g : from catalog" << endl;
  Stream << "#r : from catalog" << endl;
  Stream << "#i : from catalog" << endl;
  Stream << "#z : from catalog" << endl;
  Stream << "#ue : from catalog" << endl;
  Stream << "#ge : from catalog" << endl;
  Stream << "#re : from catalog" << endl;
  Stream << "#ie : from catalog" << endl;
  Stream << "#ze : from catalog" << endl;
  Stream << "#img : image number" << endl;
  Stream << "#star : star number in the catalog (first =1)" << endl;
  Stream << "#chi2v : chi2 of this vignette per dof " << endl;
  Stream << "#chi2pdf : chi2 of PSF photometry per dof" << endl;
  Stream << "#satur : 1 if some pixels are saturated" << endl;
  Stream << "#nsatur : number of pixels  saturated" << endl;
  Stream << "#name : dbimage name" << endl;
  Stream << "#end" <<endl;
  Stream << setprecision(12);
}

#include "calibratedstar.h"

void SimPhotFit::WriteTupleEntries(ostream &Stream, const CalibratedStar &CStar) const 
{
  const string band = RefImage().Band();
  double mag=0, mage=0;    
  if (band=="u") {mag = CStar.u; mage=CStar.ue;};
  if (band=="g") {mag = CStar.g; mage=CStar.ge;};
  if (band=="r") {mag = CStar.r; mage=CStar.re;};
  if (band=="i") {mag = CStar.i; mage=CStar.ie;};
  if (band=="z") {mag = CStar.z; mage=CStar.ze;};

  int img_count = 0;

  for (VignetteCIterator i = vignetteList.begin(); i != vignetteList.end();++i, img_count++)
    {
      const Vignette *v = *i;

      string dbim_name= v->Name();

      double sigPosX = 0, sigPosY = 0;
      int posIndex = PosIndex();
      if (posIndex >= 0)
	{
	  sigPosX = sqrt(A(posIndex,posIndex));
	  sigPosY = sqrt(A(posIndex+1,posIndex+1));
	}

      int fluxIndex = FluxIndex(v);
      if (fluxIndex < 0) continue;
      double sigFlux=sqrt(A(fluxIndex,fluxIndex));

      double sigSky = 0;
      int skyIndex = SkyIndex(v);
      if (skyIndex >= 0) sigSky = sqrt(A(skyIndex,skyIndex));

      int nterms = v->NTerms();
      double chi2Vignette = -1;
      if (nterms) chi2Vignette = v->Chi2()/double(nterms);
      
      double chi2Glob;
      int ndof;
      CumulateChi2(chi2Glob, ndof, false);
      chi2Glob = chi2Glob/ndof;

      Point fittedPos( ObjectPosInImage());
      Stream << fittedPos.x << ' ' << fittedPos.y << ' '
	     << sigPosX  << ' ' << sigPosY << ' '
	     << v->GetFlux() << ' ' << sigFlux << ' '
	     << v->GetSky() << ' ' << sigSky << ' ' 
	     << v->MJD() << ' ' 
	     << v->Seeing() << ' '
  	     << v->ExpTime() << ' '
       	     << v->PhotomRatio() << ' ' 
	     << v->GFSeeing() << ' '
  	     << v->SESky() << ' '
       	     << v->SIGSky() << ' '
	     << mag << ' ' << mage << ' '
	     << CStar.ra << ' ' << CStar.dec << ' '
	     << CStar.x << ' ' << CStar.y << ' '
	//mags
	     << CStar.u << ' ' << CStar.g << ' ' << CStar.r << ' ' 
	     << CStar.i << ' ' << CStar.z << ' '
	// mags uncertainties
	     << CStar.ue << ' ' << CStar.ge << ' ' << CStar.re << ' ' 
	     << CStar.ie << ' ' << CStar.ze << ' '
	     << img_count << ' '
	     << CStar.id << ' ' 
	     << chi2Vignette << ' ' << chi2Glob << ' '
	     << v->has_saturated_pixels << ' ' <<  v->n_saturated_pixels << ' '
	     << ' '<< dbim_name 
	     << endl;
    }
}

Mat SimPhotFit::FillNightMat() const {

  
  double mjd; // modifiedjulian day
  vector<double> nightdates;
  bool isinnight;
  double timediff = 10./24.; // 10 hours
  int nimages = 0;
  Mat nightMat;
  
  for (VignetteCIterator itVig = vignetteList.begin(); itVig != vignetteList.end(); ++itVig) 
    {
      const Vignette *v = *itVig;
      if (fluxMap.find(v) != fluxMap.end())
	{ 
	  if((toDo&FIT_FLUX)&&((*itVig)->mightFitFlux)) 
	    {
	    mjd = (*itVig)->ri.ModifiedJulianDate();
	    nimages++;
	    isinnight=false;
	    for(unsigned int day=0;day<nightdates.size(); ++day) 
	      {
	      if(fabs(mjd-nightdates[day])<timediff) 
		{
		isinnight = true;
		break;
		}
	      }
	    if(!isinnight) 
	      {
		nightdates.push_back(mjd);
	      }
	    }
	}
    } 
  int nnights = nightdates.size(); // number of nights for this lightcurve
  nightMat.allocate(nnights,nimages); // szie of matrix
  // now we fill this matrix
    //        nights
  //      <----------->
  // i   |   1 
  // m   |   1
  // a   |   1
  // g   |    1
  // e   |    1
  // s   |     1 ...  
  
  int im = 0;
  for (VignetteCIterator itVig = vignetteList.begin(); itVig != vignetteList.end(); ++itVig)
    {
       const Vignette *v = *itVig;
       if (fluxMap.find(v) != fluxMap.end())
	{
	  if((toDo&FIT_FLUX)&&((*itVig)->mightFitFlux))
	    {
	      mjd = (*itVig)->ri.ModifiedJulianDate();
	      for(unsigned int day=0;day<nightdates.size(); ++day)
		{
		  if(fabs(mjd-nightdates[day])<timediff)
		    {
		      nightMat(day,im)=1;
		      break;
		    }
		}
	      im++;
	    }
	}
    }
  //cout << "========= NightMat ==========" << endl;
  //cout << NightMat << endl;
  return nightMat;
}
