#include <iostream>
#include "senearstar.h"
#include "transformedimage.h"
#include "transformedimageaddfakes.h"
#include "fileutils.h"
#include "fitsimage.h"

//img1 devra etre l'image avec le plus mauvais sing
// On doit donner en entree le nombre de Sn a ajouter


//!Attention il faut donner une SENearStarList pleine en entree
TransformedImageAddFakes::TransformedImageAddFakes(const string &Name, 
		        const ReducedImage *source, 
			const ImageTransfo *transfo, 
			SENearStarList *NSList, 
	   		const Gtransfo *oneto2) 
  : TransformedImage(Name, *source, transfo)
{
  NearStarRef = NSList;
  One2Two = oneto2->Clone();
}


TransformedImageAddFakes::~TransformedImageAddFakes() 
{ 
  delete One2Two; // delete NearStarRef;
}


ReducedImage *TransformedImageAddFakes::Clone() const
{
  return new TransformedImageAddFakes(Name(), Source(), Transfo(), NearStarRef, One2Two);
} 

bool TransformedImageAddFakes::MakeFits() 
{
  
  string fileName = FitsName();
  if (FileExists(fileName)) return true;
  cout << " Transforming (+add fakes) image "<< Name() << endl;
  FitsImage inFits(Source()->FitsName());
  
  //Ajout des Sne aux Fits avant la transformation   
  cout << " on va ajouter " <<  NearStarRef->size() << " SN a " << Source()->Name() << endl;
  NearStarRef->AddToFitsImage(inFits,One2Two);
  double defaultVal = Source()->BackLevel();
  FitsImage outFits(fileName,FitsHeader(inFits));
  Transfo()->TransformImage(inFits, outFits, Source(), this, defaultVal);
  cout << " TransformedImageAddFakes is produced. Updating parameters. " << endl;
  //   Transfo()->Update(*this);
  return true;
}
