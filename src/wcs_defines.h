#define _wcs_h_
#define WCS_PIX -1	/* Pixel WCS */
#define WCS_LIN  0	/* Linear projection */
#define WCS_AZP  1	/* Zenithal/Azimuthal Perspective */
#define WCS_TAN  2	/* Gnomonic = Tangent Plane */
#define WCS_SIN  3	/* Orthographic/synthesis */
#define WCS_STG  4	/* Stereographic */
#define WCS_ARC  5	/* Zenithal/azimuthal equidistant */
#define WCS_ZPN  6	/* Zenithal/azimuthal PolyNomial */
#define WCS_ZEA  7	/* Zenithal/azimuthal Equal Area */
#define WCS_AIR  8	/* Airy */
#define WCS_CYP  9	/* CYlindrical Perspective */
#define WCS_CAR 10	/* Cartesian */
#define WCS_MER 11	/* Mercator */
#define WCS_CEA 12	/* Cylindrical Equal Area */
#define WCS_CPS 13	/* Conic PerSpective (COP) */
#define WCS_COD 14	/* COnic equiDistant */
#define WCS_COE 15	/* COnic Equal area */
#define WCS_COO 16	/* COnic Orthomorphic */
#define WCS_BON 17	/* Bonne */
#define WCS_PCO 18	/* Polyconic */
#define WCS_GLS 19	/* Sanson-Flamsteed (GLobal Sinusoidal) */
#define WCS_PAR 20	/* Parabolic */
#define WCS_AIT 21	/* Hammer-Aitoff */
#define WCS_MOL 22	/* Mollweide */
#define WCS_CSC 23	/* COBE quadrilateralized Spherical Cube */
#define WCS_QSC 24	/* Quadrilateralized Spherical Cube */
#define WCS_TSC 25	/* Tangential Spherical Cube */
#define WCS_NCP 26	/* Special case of SIN */
#define WCS_DSS 27	/* Digitized Sky Survey plate solution */
#define WCS_PLT 28	/* Plate fit polynomials (SAO) */
#define WCS_TNX 29	/* Gnomonic = Tangent Plane (NOAO with corrections) */
#define WCS_GALACTIC	3	/* Galactic longitude and latitude */
#define WCS_ECLIPTIC	4	/* Ecliptic longitude and latitude */
#define WCS_ALTAZ	5	/* Azimuth and altitude/elevation */
#define WCS_LINEAR	6	/* Linear with optional units */
#define WCS_NPOLE	7	/* Longitude and north polar angle */
#define WCS_SPA		8	/* Longitude and south polar angle */
#ifndef PI
#define PI	3.141592653589793238462643
#endif
#define degrad(x)	((x)*PI/180.)
#define raddeg(x)	((x)*180./PI)
#define hrdeg(x)	((x)*15.)
#define deghr(x)	((x)/15.)
#define hrrad(x)	degrad(hrdeg(x))
#define radhr(x)	deghr(raddeg(x))
#define  TNX_CHEBYSHEV    1
#define  TNX_LEGENDRE     2
#define  TNX_POLYNOMIAL   3
#define	TNX_XNONE	0	/* no x-terms (old no) */
#define	TNX_XFULL	1	/* full x-terms (new yes) */
#define	TNX_XHALF	2	/* half x-terms (new) */
