#include <iostream>
#include <fstream>
#include <iomanip>

#include <fileutils.h>
#include <dicstar.h>
#include <fitsimage.h>
//#include <basestar.h>
//#include <gaussianfit.h>
#include <polokaexception.h>

static void usage(const char *pgname)
{
  cerr << pgname << " -l <zpstarlist> -s <scatterflat> -g <grid> -r <referenceimage> ( -o <outputzpstarlist> )" << endl;
  cerr << " or :  -x <x> -y <y> -s <scatterflat> -g <grid> -r <referenceimage>" << endl;
  exit(1);
}


int main(int argc, char **argv)
{
  
  string zpstarlist_filename = "";
  string scatterflat_filename = "";
  string grid_filename = "";
  string referenceimage_filename = "";
  string output_zpstarlist_filename = "grid_corrected_zpstars.list";
  float x_coord = -1000;
  float y_coord = -1000;
  
  if (argc < 9)  {usage(argv[0]);}

  for (int i=1; i<argc; ++i)
    {
      char *arg = argv[i];
      if (arg[0] != '-')
	{
	  cerr << "unexpected parameter " << arg << endl;
	  usage(argv[0]);
	}
      switch (arg[1])
	{
	case 'x' : x_coord = atof(argv[++i]); break;
	case 'y' : y_coord = atof(argv[++i]); break;
	case 'l' : zpstarlist_filename = argv[++i]; break;
	case 's' : scatterflat_filename = argv[++i]; break;
	case 'g' : grid_filename = argv[++i]; break;
	case 'r' : referenceimage_filename = argv[++i]; break;
	case 'o' : output_zpstarlist_filename = argv[++i]; break;
	default : 
	  cerr << "unknown option " << arg << endl;
	  usage(argv[0]);
	}
    }
  
  if( ( (zpstarlist_filename == "") && ( (x_coord<0) || (y_coord<0) ) )
      || (scatterflat_filename == "")
      || (grid_filename == "")
      || (referenceimage_filename == "") )
    {usage(argv[0]);}
  
  try {
    if(x_coord<0) // use zpstarlist
      if ( ! FileExists(zpstarlist_filename) )
	throw(PolokaException("cannot read "+zpstarlist_filename));
    if ( ! FileExists(scatterflat_filename) )
      throw(PolokaException("cannot read "+scatterflat_filename));
    if ( ! FileExists(grid_filename) )
      throw(PolokaException("cannot read "+grid_filename));
    if ( ! FileExists(referenceimage_filename) )
      throw(PolokaException("cannot read "+referenceimage_filename));
  
    FitsImage scatter(scatterflat_filename);
    FitsImage grid(grid_filename);
  
#define STD_NX   2048
#define STD_NY   4612


  int x_offset,y_offset;
  {
    FitsImage referenceimage(referenceimage_filename);
    x_offset = (referenceimage.Nx()-STD_NX)/2;
    y_offset = (referenceimage.Ny()-STD_NY)/2;
  }
  int scatter_nx = scatter.Nx();
  int scatter_ny = scatter.Ny();
  int grid_nx = grid.Nx();
  int grid_ny = grid.Ny();
  
  int scatter_dx = STD_NX/scatter_nx;
  int scatter_dy = STD_NY/scatter_ny;
  int grid_dx = STD_NX/grid_nx;
  int grid_dy = STD_NY/grid_ny;
  
  
  cout << "gridcorrection DEBUG offsets = " << x_offset << " " << y_offset << endl;

  if(x_coord>=0 && y_coord>=0) {
    int x_scatter = min(max(0,int( floor((x_coord-x_offset)/scatter_dx) )),scatter_nx-1);
    int y_scatter = min(max(0,int( floor((y_coord-y_offset)/scatter_dy) )),scatter_ny-1);
    int x_grid = min(max(0,int( floor((x_coord-x_offset)/grid_dx) )),grid_nx-1);
    int y_grid = min(max(0,int( floor((y_coord-y_offset)/grid_dy) )),grid_ny-1);
    
    float scatter_value = scatter(x_scatter,y_scatter);
    float grid_value = grid(x_grid,y_grid);
    float corr = scatter_value/grid_value;

    cout << "CORRECTION " << corr << endl;
    return EXIT_SUCCESS;
  }



  // now read/write zpstars
  // ---------------------
  DicStarList zpstarlist(zpstarlist_filename);
  ofstream os(output_zpstarlist_filename.c_str());
  os << setprecision(10);
  int n_stars = zpstarlist.size();
  if ( n_stars == 0 )
     throw(PolokaException(zpstarlist_filename+" is empty"));
  
  int n=0;
  bool header_written = false;
  
  for(DicStarIterator entry=zpstarlist.begin();entry!=zpstarlist.end();++entry,++n) {

    if ( (*entry)->getval("nmeas") < 10 ) {
      cerr << "only " << (*entry)->getval("nmeas") << " measurements for this star, skip it" << endl;
      continue;
    }
    
    

    float x  = (*entry)->getval("x");
    float y  = (*entry)->getval("y");
    
    if(x==0 || y==0) {
      // need to be fix with exception in dicstar
      throw(PolokaException("I find x or y = 0, wrong key in file??"));
    }

    float flux = (*entry)->getval("flux");
    
    int x_scatter = min(max(0,int( floor((x-x_offset)/scatter_dx) )),scatter_nx-1);
    int y_scatter = min(max(0,int( floor((y-y_offset)/scatter_dy) )),scatter_ny-1);
    int x_grid = min(max(0,int( floor((x-x_offset)/grid_dx) )),grid_nx-1);
    int y_grid = min(max(0,int( floor((y-y_offset)/grid_dy) )),grid_ny-1);
    
    float scatter_value = scatter(x_scatter,y_scatter);
    float grid_value = grid(x_grid,y_grid);
    float corr = scatter_value/grid_value;
    float corrected_flux = flux*corr;
    (*entry)->setval("flux",corrected_flux);
    cout << x << " " << y << " " << x_scatter << " " << y_scatter << " " << corr << endl;
    
    (*entry)->AddKey("gridcorr",corr);

    if(! header_written) {
      (*entry)->WriteHeader_(os);
      header_written=true;
      os << "#end" << endl;
    }
    (*entry)->write(os);
    
    

  }
  os.close();

  // a developper



  }catch(PolokaException p)
    {
      p.PrintMessage(cout);
      return EXIT_FAILURE;	      
    }
  
  return EXIT_SUCCESS;
}
