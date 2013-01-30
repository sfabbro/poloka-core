#include "poly1.h"

#include <fstream>

using namespace std;

void Poly1::Write(ostream &file) const
{
  file << " Poly1 1 " << endl;
  file << npar << endl;
  file << a << ' ' << b << endl;
  for (unsigned k=0; k<npar; ++k) file << params(k) << ' ';
  file << endl;
}


void Poly1::Write(const string &FileName) const
{
  ofstream file(FileName.c_str());
  Write(file);
  file.close();
}


bool Poly1::Read(ifstream &file, const string &FileName)
{
  string type;
  int version;
  file >> type >> version;
  if (type != "Poly1" || version != 1)
    {
      cout << " wrong file type and or version in file " << FileName << endl;
      file.close();
      return false;
    }
  file >> npar;
  file >> a >> b;
  params.allocate(npar);
  for (unsigned k=0; k<npar; ++k) file >> params(k);  
  return true;
}



bool Poly1::Read(const string &FileName)
{
  ifstream file(FileName.c_str());
  if (!file)
    {
      cout << " cannot open " << FileName << endl;
      return false;
    }
  bool status = Read(file, FileName);
  file.close();
  return status;
}

void Poly1::operator += (const Poly1 &Right)
{
  if (npar != Right.npar) 
    {
      std::cout << " Poly1::operator += : different degrees " << std::endl;
      abort();
    }
  params += Right.params;
}

