#ifndef SCORESCOLLECTION__H
#define SCORESCOLLECTION__H

#include <string>
#include <vector>
#include <iostream>

class ScoresCollection;

ScoresCollection* FindCollection(const std::string &C);


void StoreScore(const std::string &Context, const std::string &Name, 
		const std::vector<double> & Vals);

void StoreScore(const std::string &Context, const std::string &Name, 
		const double Value);

void WriteScores(const std::string &Context, std::ostream &S, 
		 const bool WithContext=false);

void WriteScores(const std::string &Context, const std::string &FileName, 
		 const bool WithContext=false);

void WriteAllScores(std::ostream &S, const bool WithContext=false);

void WriteAllScores(const std::string &FileName, const bool WithContext=false,
		    const bool Append = false);


#endif
