#ifndef RESAMPLER__H
#define RESAMPLER__H

/* This file implements what concerns the resampling for the simultaneous
photometric fit (the new version without resampling prior to fitting).
The basic image resampling routines are not exactly what is nedded here,
since we need derivatives of the output image w.r.t the input image,
which is not needed when just resampling an image.

   We need basicaly 2 or 3 funtionalities:
- compute the derivatives
- get the "size overhead" (the number of lost pixel on image boundaries)
- not mandatory : image resampler that matches the derivatives. We may
  implement it for sake of efficiency.

These 2 or 3 functionalities have to be consistent.
If you consider changing something somewhere, check that the consistency
is preserved (for example, that the derivatives of the provided
resampler are the ones you claim.
*/

class Gtransfo;
class Array4D;

void ResamplerComputeDerivatives(const Gtransfo *Transfo, 
				 Array4D & Coeffs);

int ResamplerBoundarySize();

class PixelBlock;
//! does not assign boundaries of "Out"
void ResampleImage(const PixelBlock &In, const Gtransfo* Out2In,
		   PixelBlock &Out);

//! Assigns boundaries of Out
void ConvolveImage(const PixelBlock &In, const PixelBlock &Kernel,
		   PixelBlock &Out);

#endif /* RESAMPLER__H */
