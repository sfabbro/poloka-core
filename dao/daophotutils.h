// This may look like C code, but it is really -*- C++ -*-
#ifndef DAOPHOTUTILS__H
#define DAOPHOTUTILS__H

class ReducedImage;
class Daophot;
class PsfStars;

//! Iterate building of a daophot PSF 
bool IteratePsf(Daophot& DaoSession, PsfStars& Stars);

//! Produce a PSF using DAOPHOT
void MakeDaoPsf(ReducedImage &Rim, const bool Redo=false);

void MakeExperimentalPsf(ReducedImage &Rim);

//! Produce a PSF and a ALLSTAR catalog. Can merge the result catalog into an existing se.list
void MakeDaoPsfAls(ReducedImage &Rim, const bool Merge, const bool Redo=false);

//! Build a fake stars using DAOPHOT routines and PSF
void MakeFakeStarImage(const ReducedImage &Rim, const bool Redo=false);

//! Update a ReducedImage with its associated DAOPHOT PSF file
bool UpdateSeeingFromDaoPsf(ReducedImage &Rim);


#endif // DAOPHOTUTILS__H
