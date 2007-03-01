c     MAXBOX is the square subarray that will hold the largest final PSF.
c     If the maximum PSF radius permitted is R, then MAXBOX is the
c     odd integer 2*INT(R)+1.  However, because we will be dealing
c     with two levels of interpolation:  (1) interpolating the raw
c     picture data to arrive at a PSF whose centroid coincides with
c     the central pixel of the look-up table; and (2) interpolating
c     within the PSF itself to evaluate it for comparison with the
c     raw picture data for the program stars, PSF will have to
c     operate on a square array which is larger by 7 pixels in X
c     and Y than MAXBOX.  Hence, the dimensions below are all
c     MAXBOX + 7.
      INTEGER MAXBOX
      PARAMETER(MAXBOX=69)

c     MAXPSF is the dimension of the largest lookup table that will ever 
c     need to be generated.  Recall that the corrections from the
c     Gaussian approximation of the PSF to the true PSF will be
c     stored in a table with a half-pixel grid size.
c     MAXPSF must then equal     
c     2*[ 2*INT(PSFRADMAX) + 1 ] + 7.
      INTEGER MAXPSF
      PARAMETER(MAXPSF=2*(2*MAXBOX+1)+7)

c     MAXPAR is the maximum number of parameters which may be used to 
c     describe the analytic part of the PSF, *** IN ADDITION TO
c     the central intensity, and the x and y centroids.
      INTEGER MAXPAR
      PARAMETER(MAXPAR=6)

c     MAXEXP is the maximum number of terms in the expansion of the
c     look-up table (1 for constant, 3 for linear, 6 for quadratic).
      INTEGER MAXEXP
      PARAMETER(MAXEXP=10)
      
c     MAXSTR is the largest number of stars permitted in a data file.
      INTEGER MAXSTR
      PARAMETER(MAXSTR=200)
c this number cannot be too large because there are arrays which are 
c sized as (3*MAXSTR,3*MAXSTR). Pierre Astier.


      
c     MAXN is the maximum permitted number of PSF stars.
      INTEGER MAXN
      PARAMETER(MAXN=200)
      
c     MAXSKY is the maximum number of sky pixels we can deal with,
c     given the limited amount of working space.      
      INTEGER MAXSKY
      PARAMETER(MAXSKY=10000)

c     MAXMAX is the largest group for which a solution will ever be
c     attempted = maximum permissible value of MAXGRP.     
      INTEGER MAXMAX
      PARAMETER(MAXMAX=100)
      
      INTEGER NOPT
      PARAMETER(NOPT=20)
      
  
