// We use the casio graphics library,2018
// documentation: https://www.cairographics.org/manual/cairo-Paths.html
#include <unistd.h>
#include <assert.h>
#include <cairo/cairo.h>
#include <cairo/cairo-xlib.h>
#include <X11/Xlib.h>
#include <X11/Xutil.h>
#include <stdio.h>

typedef struct {
  double x,y,vx,vy,Rad,Mass,Contamination_time;
  int State; // 0 -> particule saine, 1 -> particule infectée, 2 -> particule rétablie, 3 -> mur
  int Zone;
} Particle; //this is our definition of  a particle, position plus speed


class Graphics{
 private:
  int Np;
  double lmin_x, lmin_y, lmax_x, lmax_y,DimX,DimY,diam, Ratio;
  cairo_surface_t *sfc;
  cairo_t *cr;
  Display *dsp;
  Drawable da;
  int count;
 public:
  Graphics(int N, int Pix,double Lmin_x,double Lmin_y, double Lmax_x,double Lmax_y, double diameter){
    Np=N;
    count=1000;
    assert( Lmin_x < Lmax_x );
    assert( Lmin_y < Lmax_y );
    assert(Np >0) ;
    lmin_x=Lmin_x;
    lmin_y=Lmin_y;
    lmax_y=Lmax_y;
    lmax_x=Lmax_x;
    Ratio = (Lmax_x - Lmin_x) / ( Lmax_y - Lmin_y );
    double PixArea = Pix*Pix;
    DimY=sqrt(PixArea/Ratio);
    DimX=Ratio*DimY;
    diam=diameter;

    if (( dsp = XOpenDisplay(NULL)) == NULL)      exit(1);//window management X11, and cairo graphics
    int screen = DefaultScreen(dsp);
    da = XCreateSimpleWindow(dsp, DefaultRootWindow(dsp), 0, 0, DimX, DimY, 0, 0, 0);
    XMapWindow(dsp, da);
    sfc = cairo_xlib_surface_create(dsp, da, DefaultVisual(dsp, screen), DimX , DimY);
    cairo_xlib_surface_set_size(sfc, DimX, DimY);
    cr = cairo_create (sfc);
  } //empty window now on screen

  //some error messages on trying to duplicate the window
  Graphics & operator=(Graphics &g){// stop the program when copying the window
    printf("Don't use = with graphics objects  %p \n", &g ); exit(1);
  }
  Graphics(const Graphics &g ){// stop the program when passing windows as argument
    printf("Don't pass graphics objects: use a pointer  %p \n" , &g); exit(2);
  }
  void writePNG(){//save to png file, using "count" to produce a numbered file;
    char c[100];
    snprintf(c,99,"md%d.png", count++ );
    cairo_surface_write_to_png(sfc, c );
  }



  void draw(Particle*, double, int, double **, int**, double**, int**, double *);//draw the particles
  void frame(double , double , double , double );//draw a square
  ~Graphics(){cairo_destroy (cr);cairo_surface_destroy (sfc); } //clean up function
};
