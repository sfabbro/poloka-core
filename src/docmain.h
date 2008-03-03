/* This file is NOT to be included by any part of the code.
   It is parsed by doxygen to produce the documentation. */

/*! \mainpage The Poloka (aka TOADS) Software
        TOols for Analysis and Detection of Supernovae


   \authors Pierre Astier <BR>
            Sebastien Fabbro <BR>
            Julien Guy <BR>
            Delphine Hardin <BR>
            Julien Raux  <BR>
            Nicolas Regnault <BR>
            Kyan Schahmaneche <BR>

	       All these people belong(ed) to the <A HREF="http://supernovae.in2p3.fr/group.html">FROGS </a> <BR>



   \date \today


   \attention This document is the user manual and reference manual of Poloka. 
   It will then never be up to date. It is generated using "doxygen".<BR>
   It is far from being complete.


The structure of the document is:

<ul>
 <li> \ref introduction
 <li> \ref images
 <li> \ref imageio
 <li> \ref geomtransfo
 <li> \ref starlists
 <li> \ref database
 <li> \ref reducedimage
 <li> \ref fitstoad
 <li> \ref env_var
 <li> \ref gettingstarted
   <ul>
   <li> \ref installing_images 
   <li> \ref inspecting_database 
   <li> \ref cataloging
   </ul>
 <li> \ref subtracting
 <li> \ref lightcurve
</ul>

\section introduction Introduction

TOADS aims at providing a framework for astronomy image analysis, with an
emphasis on image subtraction and differential photometry of varying objects.
At variance with most of the codes developped in the astronomy community,
it is NOT intended to provide ready-to-use tools. There is in fact
very little chance that your problem is already solved in Poloka,
unless you are doing the same thing as we are. Our purpose
is more to provide a convenient toolbox, in which you
will have to collect components and write your own code.
If C++ evocates nothing to you, Poloka is likely
to be the wrong choice.

   Our second bias is that there is nothing interactive (yet?) in TOADS,
because it was developped originally to carry out supernova
detections by image subtraction, while the observers are asleep.
And it actually does that, and more things as time goes on.

  We have tried to separate as much as possible the various concepts
involved in image analysis, such as images, geometrical transformations,
lists of stars.

\section images Images
the images are implemented in memory as arrays of float, with 
pixel addressing (starting at 0) and most of the algebra 
operators are implemented. see Image.

\section imageio Images I/O
The Image 's are read/written using cfitsio. The access to header keys
are in the FitsHeader class, the FitsImage class assembles a FitsHeader
and an Image. read and write (if applicable) are associated with
constructors/destructor. There is no check in the provided code 
that you are altering the contents associated with a file opened read-only,
i.e. that you are modifying data in memeory and only in memory, without
any effect on disk.


  There are 2 auxialiary classes for sets of fits images: FitsSet,
which checks that the images where taken with the same filters and instrument,
and corresponds to the same cgip, and FitsParallelSlices, which is intended to 
allow processing a large number of images in a very parallel way,
ass needed for e.g. superflat computation or image stacking.

\section geomtransfo Geometrical Transformations


  The concept of geometrical transformation mainly refers
to applications that transform a point of the plane into a point
of the plane. We then have routines that transform images and 
lists of stars.

  The abstract concept was implemented as an abstract class : Gtransfo.
The actual implementations that already exist are polynomial
transformations of inscreasing degree: GtransfoLin, GtransfoQuad, GtransfoCub,
rotations and translations are defined from GtransfoLin, and 
the do-nothing GtransfoIdentity, which is handy in several cases.
  On top of providing the transformation, the classes should provide
a fitting routine, which uses a collection of pairs of points to
be fitted to the model the class implements. These list are called
StarMatchList and can be constructed in various ways. Directly
from starlists (see below), you can invoque the guessing routines
in listmatch.h. There are also higher level wrapping routines 
in imagematch.h


 \section starlists Stars and lists of stars
A few star formats are already defined in the code. The BaseStar only contains 
coordinates and flux, and is mainly used for geometrical matching. SEStar is the 
star that we extract from SEXtractor, and it derives from BaseStar. All star kinds
have a StarList defined, which uses the STL list container. They are coded
in such a way that lists of derived objects can be casted to lists of a base type.
This is what the SE2Base() routines do.


 \section database Database
 We have implemented a very crude file database, which is in fact not 
a database. An entry is this database, is a DbImage, which is an image in an abstract sense:
a DbImage contains a raw image, a flat image, a dead image, a fringe image, 
a catalog, and other goodies. What the Database does is mainly to provide 
where the datafiles are, provided a name (the DbImage name) and a function.
The database code searches the requested information in directories 
which are to be provided by the user. See \ref database_page for details.


 \section reducedimage Various image types
  The successive steps for Image subtraction are handled through a set
of classes derived from ReducedImage: TransformedImage, ImageSum, ImageSubtraction.
They all follow the same pattern as DbImage (they in fact inherit from it).
When you construct such images (you transform, sum, subtract), ReducedImage's
are created and written to disk. The physical location is the
one tagged by "here" in the \red dbconfig. The actual value is usually
the local directory (As in \ref dbconfig_example), from which you actually run
your executable (e.g. newmake_sub).


\section gettingstarted Getting Started.

To obtain the code, you have to go through one of the authors.
Installing the code and prerequisite is detailed in the README.
before running anything that uses DbImage's, you have to provide a db config file
(see \ref dbconfig). make_catalog and image subtractions use datacards,
which are located with TOADCARDS, which has to be defined as the directory
that contains those datacards (distributed with the source code).

\subsection installing_images Putting images in the data base.
  You have to use install_image for that.
  If you want to install the image 604053o02.fits in the directory
  /snovad24/cfht01B, you just have to do:
\code
 install_image /snovad24/cfht01 604053o02.fits
\endcode
You can install a bunch of images at a time. The DbImage
you just created will be named 604053o02. It now has one
field: the "raw" image. If the fits image you want to enter
is already flatfielded, you should invoke install_image
  with the option "-cal". In this case, the assigned field of your
DbImage is the "calibrated" image. Notice that the fits file
is not copied but linked. The routine that assigns this link does
  it best to avoid absolute pathes in the links, in order to allow
you to move a whole directory tree if needed.

  Then for the database to find the image back, you have to put
the path /snovad24/cfht01 in your \ref dbconfig.
  When you are done, you can check that the image is actually located
  by the database system:

\subsection  inspecting_database Inspecting Database.
  The command dbls enables you to see what the database currently contains.
  dbls -h issues a small help. dbls is far less sophisticated than ls.
  "dbls 604053o02" should output 604053o02, and "dbls -raw 604053o02"
should output /snovad24/cfht01/604053o02/raw.fits. This is the file name
that you will get in your code when you will want to access the "raw" field
of your DbImage.
  We have a "header" command that enables you to dump a fits header or 
subparts of it, and it is connected to the database:
\code
> header -k RA -k DEC /snovad24/cfht01B/604053o02/calibrated.fits
/snovad24/cfht01B/604053o02/calibrated.fits 22:17:50.59 0:10:26.3
> header -k RA -k DEC -cal 604053o02
604053o02 22:17:50.59 0:10:26.3 
\endcode

\section fitstoad Standardization of headers.
Sebastien Fabbro has designed a scheme that enables to present the fits headers
in a telescope/instrument independent fashion to the code. This is based 
  on active routines, rather than overwritting actual headers. This goes
  through a virtual instrument system, that performs the right translations
and/or computations for a specific telescope instrument. There is a default
implementation that may do correctly for a telescope/instrument still 
unknown to TOADS. You may test it through:
\code
header -k TOADALL your_image.fits
\endcode
To see how well (or badly) it works for your specific case. If it does 
not work properly, you have to go through \ref newtelinst. Drop us
a mail if you get lost.


\section env_var Environment variables.
  Various parts of the code rely (or may depend) on envirnment variables.
  Here is a tentative list of those:
  <ul>
  <li> DBCONFIG : where the "dbconfig" (see \ref dbconfig) file stands. Default
         location is \verb "~/.dbconfig". Providing such a file is mandatory
       to use any code that relies on more than a single fits file. see \ref dbconfig_example for and example.
       <li>  TOADSCARDS : the directory where the datacards are 
        (somewhere/poloka/datacards). This directory contains 
       the SExtractor usual configuration files.  The Poloka tunable 
       quantities are in sub.datacards (with hopefully meaningful comments...).
       <li> USNODIR, USNOFILE. Poloka provides a functionality to match 
  an image to a reference catalog (matchusno). This reference catalog 
   is the USNO by default ands its location  is then provided through USNODIR.
   You may as well provide your own catalog on the command line or through
   USNOFILE.  
  </ul>

\subsection cataloging Making a catalog.
Most of the useful code of Poloka requires image catalogs. You just
have to run make_catalog to get it computed and written to disk.
  This constructs a SExtractor catalog, and many other goodies
  such as a weight map (from sky variance and flat if available),
  a saturated pixel map, a tentative detection of cosmics and
  satellite trails (incorporated to the weight map).


\subsection matchusno Computing a WCS
Poloka contains code to match an image to a reference catalog.
  It actually matches the catalog of an image to a reference catalog,
  and deduces a WCS, and writes it into the image fits header.
The starting point of the match relies on the user providing
  a guess, through a WCS. This starting WCS can be either 
  provided in the image header (most of the observing systems now
  compute and write a decent WCS from telescope settings), or computed
  using keywords of the header. This is done on a telescope-instrument basis
  in the "telinst" directory (see \ref newtelinst) by providing routines
  that will compute a rough WCS for every image (without actually 
  writing it anywhere). The guess WCS is mainly used to collect the
  right area in the reference catalogue. Then a shift is tried,
  and when it fails, rotations (with and without flip) are attempted. 
  To see how well or bad the code is doing, "matchusno -n " does everything
  but writing the found WCS. Our WCS follow the standards, and encode
  the optical distortions in a way Swarp understands (which has perhaps
  become a standard by now). The reference catalog can be provided
  as the USNO one (2.0), or as an external catalog (see \ref env_var).
  An external catalog can be provided on the matchusno command line
  (see \ref  usno_file_format if you need)


\section subtracting Running a subtraction.
You have to create a "subfile" (see \ref subfile), 
and run makesub. Input images should not be aligned in advance,
but they should overlap. If you really need to mosaic images,
you should definitly contact us. 
  You may log the output, which contains very useful informations when things
go wrong.


\section lightcurve Producing a lightcurve
To build a lightcurve for a supernova we need to know where is the 
supernova and on which night it has been observed. You then need 
to produce a "lightfile" (see \ref lightfile) which mainly
  describes which images to use, and which contain light of the variable 
  objet. These images must be geometrically aligned (i.e. be on the same
  pixel grid).
Find a directory with a lot of space. Then type

     make_lightcurve <myfile>

and go for coffee. See \ref lcresults for lightcurve results.
This is basicaly the code used (on the French side) to produce the 
SNLS light curves.
A new scheme for computing light curves from  unregistered images has 
been developed, and   passed the first tests. At variance with the 
one decribed above, it was never tested on a large scale.



*/
