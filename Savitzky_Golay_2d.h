#ifndef SAVITSKY_GOLAY_2D_H
#define SAVITSKY_GOLAY_2D_H
#include <fftw3.h>
#include <Accelerate/Accelerate.h>
using namespace std;

class Range2d
{
public:
    Range2d(int xsize, int ysize){
        xsize_ = xsize;
        ysize_ = ysize;
        xsize_1 = xsize-1;
        ysize_1 = ysize-1;
    }
    Range2d(){}
    int init(int xsize, int ysize){
        xsize_ = xsize;
        ysize_ = ysize;
        xsize_1 = xsize-1;
        ysize_1 = ysize-1;
        return 0;
    }
    size_t operator()(int i, int j)
    {
        while(i < 0) i+= xsize_;
        while(i > (xsize_1)) i -= xsize_;
        while(j < 0) j+= ysize_;
        while(j > (ysize_1)) j -= ysize_;
        return i + j * xsize_ ;
    }
protected:
    int xsize_;
    int ysize_;
    int xsize_1;
    int ysize_1;
};


class Filter2d
{
public:
    Filter2d();
    Filter2d(int nx, int ny, int ld, int m, int sizex, int sizey);
    void reset(int nx, int ny, int ld, int m);
    int apply(double *indata, double *outdata, int datas_size);
    ~Filter2d(){}
protected:
    Range2d box_range;
    Range2d filter_range;
    int x_size;
    int y_size;
    int size;
    int size2;
    int roi;
    int nx;
    int ny;
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
