#ifndef Savitzky_Golay_3d_h
#define Savitzky_Golay_3d_h
#include <fftw3.h>
#include <Accelerate/Accelerate.h>
#include "range.h"
using namespace std;


class Range3d
{
public:
    Range3d(int xsize, int ysize, int zsize){
        xsize_ = xsize;
        ysize_ = ysize;
 	    zsize_ = zsize;
        xsize_1 = xsize-1;
        ysize_1 = ysize-1;
	    zsize_1 = zsize-1;
	    xylen = xsize*ysize;
    }
    Range3d(){}
    int init(int xsize, int ysize, int zsize){
        xsize_ = xsize;
        ysize_ = ysize;
        zsize_ = zsize;
        xylen = xsize*ysize;
        xsize_1 = xsize-1;
        ysize_1 = ysize-1;
        zsize_1 = zsize-1;
        return 0;
    }
    size_t operator()(int i, int j, int k)
    {
        while(i < 0) i+= xsize_;
        while(i > (xsize_1)) i -= xsize_;
        while(j < 0) j+= ysize_;
        while(j > (ysize_1)) j -= ysize_;
	    while(k < 0) k+= zsize_;
     	while(k > (zsize_1)) k -= zsize_;
        return i + j * xsize_ + k * xylen ;
    }
protected:
    int xsize_;
    int ysize_;
    int zsize_;
    int xsize_1;
    int ysize_1;
    int zsize_1;
    int xylen; 
};


class Filter3d
{
public:
    Filter3d();
    Filter3d(int nx, int ny, int nz, int ld, int m, int sizex, int sizey, int sizez);
    int apply(double *indata, double *outdata, int datas_size);
    ~Filter3d(){}
protected:
    Range3d box_range;
    Range3d filter_range;
    int x_size;
    int y_size;
    int z_size;
    int size;
    int size2;
    int roi;
    int nx;
    int ny;
    int nz;
    int diff;
    int poly_order;
    double *data;
    fftw_complex *datat;
    fftw_complex *filter;
    fftw_complex *datat_filtered;
    fftw_plan plan_data2filter;
    fftw_plan plan_data2datat;
    fftw_plan plan_filtered2data;
    int alloc();
    int free();
};

#endif
