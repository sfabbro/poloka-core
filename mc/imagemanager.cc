#include <iostream>
#include <iomanip>
#include <ctime>
#include "imagemanager.h"
#include "reducedimage.h"
#include "fitsimage.h"
 
static string IntToString(int nombre,int min_char = 0 )
{
	string res="";

	if( min_char>0)
		{
		int i;
		if(nombre != 0) i = min_char - int(log10(double(abs(nombre)))+1.0);
		else i = min_char - 1 ;
		while(i>0){res = '0' + res; i--;}
 		}
	if(nombre < 0){nombre= 0-nombre;res = '-'+res;}	
	char a = nombre%10 +'0';
	if (nombre/10 !=0) res = res + IntToString(nombre/10,0);
	res=res+a;
	return (res);
};




#include "astroutils.h" // IdentifyDeepField.


one_image::one_image(const string &image_name)
{

ReducedImage current(image_name);
FitsHeader header(current.FitsName());


// ************************** name
name = image_name;

// ************************** filtre
string temp;

	temp=current.Band();
	if (temp.size()==1)
	{
		if(temp[0]> 64 && temp[0]< 91) temp[0]=temp[0]+32;
	}
	else {cerr << "Band error" << endl; exit(1);} 
filtre = temp;

// ************************** field
	
 temp = string(header.KeyVal("OBJECT"));

	if ( temp[1]>'0' &&temp[1]<'5' )
	{
		field = (int)temp[1] - (int)'0';
	}
	else 
	  {
	    if (IdentifyDeepField(header, temp)) // try using coordinates
	      {
		int f = (int)temp[1] - (int)'0';
		if (f>0 && f<5) field = f;
	      }
	    else {cerr << "field error : " << str_field << endl; exit(1);}
	  } 
str_field  = "D" + IntToString(	field);
// *************************** ccd
ccd = header.KeyVal("EXTVER");

str_ccd = "ccd_" + IntToString(ccd,2);
		
// *************************** date



date = header.KeyVal("MJDATE");

 double MJD = date;

// calcul du time_t de ref
const double MJDREF = 52914;
const int jour = 3600*24;
time_t mytime;
tm * timeinfo;
timeinfo = gmtime(&mytime);
timeinfo->tm_isdst = 1;
timeinfo->tm_year = 2003-1900;
timeinfo->tm_mon =  10-1;
timeinfo->tm_mday = 1;
timeinfo->tm_hour  =12;
timeinfo->tm_min  = 0;
timeinfo->tm_sec = 0;
const time_t SECREF= mktime(timeinfo);


// calcul et creation de la str_date



int int_MJD = (int)MJD;
mytime = (time_t)((int_MJD-MJDREF)*jour) + SECREF;
timeinfo =  gmtime(&mytime);
 str_date = IntToString(timeinfo->tm_year+1900,4)+"-"+
		  IntToString(timeinfo->tm_mon+1,2)+"-"+
		  IntToString(timeinfo->tm_mday,2);


Done = false;

#ifdef DEBUG
cout <<  name 
<<"\t\tfiltre : " <<  filtre
 
<<"\t\tfield : " << field << "   " << str_field

<<"\t\tccd  : " << ccd  <<"   "<< str_ccd
 
<<"\t\tdate : " << date<< "   " << str_date<< endl;

#endif


};



void one_image::Dump()
{
cout <<  name 


<<"\t\tfield : " << str_field

<<"\t\tdate : " << str_date

<<"\t\tfiltre : " <<  filtre
 


<<"\t\tccd  : " << str_ccd << "\t\t"<< Done
 
<< endl;



}



Image_Manager::Image_Manager(StringList& list)
{

int total = list.size();
cout << "loading list" << endl;
int i = 0;
 for (StringIterator it = list.begin(); it != list.end();it++)
 {
	i++;
	cout << i << "\t/\t" << total << endl; 
	one_image in_image(*it);
	push_back(in_image);

 }
}



 
int Image_Manager::AddToStringList(StringList& res)
{
int i =0;
  for (ManagerIterator it = begin(); it != end();)
 {

	res.push_back(it->Name());
	++it;
	i++;
 }
 return i;
}


void Image_Manager::Dump()
{


for (ManagerIterator it = begin(); it != end();it++)  it->Dump(); 

};

ManagerIterator Image_Manager::First_with( const int in_date, const string &in_filtre, const int ccd)
{
  for (ManagerIterator it = begin() ; it != end();it++)
    if((it->int_Date() == in_date || in_date == 0) && 
       (it->Filtre()==in_filtre || in_filtre =="") &&
       (ccd ==99 || ccd == it->CCD())) return it;
 return end();
}





























