#ifndef IMAGEMANAGER__H
#define IMAGEMANAGER__H

#include <string>

#include "stringlist.h"
 
 using namespace std;
 class one_image {
 private:
 string name;
 string filtre;
 
 int field;
 string str_field;

 int ccd;
 string str_ccd;
 
 double date;
 string str_date;
 
 public:
 one_image(const string &image_name);
 ~one_image(){};
 bool Done;
 
 int int_Date(){return (int)date;};
 double Date(){return date;};
 string str_Date(){return str_date;};
 
 string Name(){return name;};
 
 int CCD(){ return ccd;};
 string str_CCD(){return str_ccd;};
 
 int Field(){return field;};
 string str_Field(){return str_field;};
 
 string Filtre(){return filtre;};
 string Sub_Path(){return (str_field+"/"+str_date+"/"+filtre+"/sub/"+str_ccd);};
 string Dbim_Path(){return (str_field+"/"+str_date+"/"+filtre+"/dbim");}; 
 void Dump();
 
 
 };
 
inline bool operator < ( one_image &oper1 , one_image &oper2)
 {
 if		(oper1.Field()  != oper2.Field()) 		return oper1.Field()  < oper2.Field();
 else if	(oper1.int_Date()   != oper2.int_Date()) 	return oper1.int_Date()   < oper2.int_Date();
 else if	(oper1.Filtre() != oper2.Filtre()) 		return oper1.Filtre() < oper2.Filtre();
 else 								return oper1.CCD()    < oper2.CCD();
 
 };
 
typedef list<one_image>::iterator ManagerIterator;
typedef list<one_image>::const_iterator ManagerCIterator;


 class Image_Manager: public list<one_image> {
 
 
 private:

 
 public:
 Image_Manager(){};
 Image_Manager(StringList& list); 
int AddToStringList(StringList& res); // return number of added names
void Dump(); 
 ManagerIterator First_with(const int in_date=0, const string &in_filtre = "", const int ccd = 99);
 
 
 
 };
 
 
  
 
#endif
