#include "nonlinmodel.h"

#include "fstream"

using namespace std;

bool NonLinModel::Write(const string &FileName) const
{
  ofstream file(FileName.c_str());
  file << " NonLinModel 1 " << endl;
  file << " maxPix " << maxPix << endl;
  Poly1::Write(file);
  for (unsigned k1=0; k1<NPar(); ++k1)
    {
      for (unsigned k2=0; k2<NPar(); ++k2)
	file << cov(k1,k2) << ' ';
      file << endl;
    }
  file.close();
  return true;
}


NonLinModel::NonLinModel(const std::string& FileName) : Poly1(-1,0,1)
{
  Read(FileName);
}

bool NonLinModel::Read(const std::string& FileName)
{
  ifstream file(FileName.c_str());
  if (!file)
    {
      // Should throw a NonExistentFileException (still to be written)
      cout << "NonLinModel::Read cannot open " << FileName << endl;
      return false;
    }
  cout << " NonLinModel::Read : untested code (psf/nonlinmodel.cc) remove this print when it works " << endl;  
  string s;
  int version;
  file >> s >> version;
  if (s != "NonLinModel" || version != 1)
    {
      cout << " NonLinModel::Read : wrong class name and/or version in " << FileName << endl;
      return false;
    }  
  file >> s;
  if (s != "maxPix" )
    {
      cout << " NonLinModel::Read : wrong format in " << FileName << endl;
      return false;
    }
  file >> maxPix;
  Poly1::Read(file, FileName);
  cov.allocate(NPar(), NPar());
  for (unsigned k1=0; k1<NPar(); ++k1)
    {
      for (unsigned k2=0; k2<NPar(); ++k2)
	file >> cov(k1,k2);
    }
  file.close();
  return true;
}
