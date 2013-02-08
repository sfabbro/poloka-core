// just a rebinning
void rebin_image(Image & image, int rebinx, int rebiny) {
  if(rebinx<=0 || rebiny<=0) {
    cerr << "ERROR rebin_image " << rebinx << " " << rebiny << endl;
  }
  // memory consumming 
  int newnx = image.Nx()/rebinx;
  if(image.Nx()%rebinx)
    newnx++;
  int newny = image.Ny()/rebiny;
  if(image.Ny()%rebiny)
    newny++;
  
  Image newimage(newnx,newny);
  Image newweight(newnx,newny);
  Pixel *p = image.begin();
  Pixel *pend = image.end();
  int nx = image.Nx();
  int inx,iny;
  int ix,iy;

  ix=iy=inx=iny=0;
  for (; p!=pend;++p) {
    newimage(inx,iny) += (*p);
    newweight(inx,iny) += 1.;
    ix++;
    if(ix>=nx) {
      ix=0;
      iy++;
      iny = iy/rebiny;
    }
    inx = ix/rebinx;
    //if(inx>= newnx || iny>= newny) {
    // cout << "error " << ix << " " << iy << " = " << inx << " " << iny << endl;
    //}
  }
  newimage/=newweight;
  image=newimage;
} 


// just an unbinning (need something clever here, for the future)
void unbin_image(Image & image, int nx, int ny, int rebinx, int rebiny) {
  if(nx<=image.Nx() || ny<=image.Ny()) {
    cerr << "ERROR rebin_image " << nx << " " << ny << endl;
  }
  
  // memory consumming 
  Image newimage(nx,ny);
  
  Pixel *p = newimage.begin();
  Pixel *pend = newimage.end();
  int inx,iny;
  int ix,iy;
  
  ix=iy=inx=iny=0;
  for (; p!=pend;++p) {
    (*p) = image(inx,iny);
    ix++;
    if(ix>=nx) {
      ix=0;
      iy++;
      iny = iy/rebiny;
    }
    inx = ix/rebinx;
  }
  image=newimage;
} 

