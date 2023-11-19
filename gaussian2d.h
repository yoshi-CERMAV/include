#ifndef GAUSSIAN_2D
#define GAUSSIAN_2D
#include "mpfit.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <iostream>
#include <sstream>
using namespace std;
const double gauss_width_cof = 2.7725887222397811449070559319807216525077819824219;
const double gauss_area_cof = 0.93943727869965132359908466241904534399509429931641;  //sqrt(ln2/pi)*2

class Gaussian2D
{
public:
    Gaussian2D(){}
    ~Gaussian2D(){}
    void set(const double *p){
        height = p[0]; pos_x = p[1]; pos_y = p[2]; width_x = p[3]; width_y = p[4];
    }
    inline double jacobian(const double &x, const double &y, double **der, int i)
    {
        double arg_x = (x-pos_x)/width_x;
        double arg2_x = arg_x*arg_x;
        double arg_y = (y-pos_y)/width_y;
        double arg2_y = arg_y*arg_y;
        double arg2 = arg2_x+arg2_y;
        double val = exp(-gauss_width_cof*arg2);
        double fac_x = 2 * arg_x / width_x;
        double fac_y = 2 * arg_y / width_y;
        if(der[0]) der[0][i] = val;
        double dpx = fac_x * gauss_width_cof * val;
        double dpy = fac_y * gauss_width_cof * val;
        if(der[1]) der[1][i] = dpx;
        if(der[2]) der[2][i] = dpy;
        if(der[3]) der[3][i] = dpx * arg_x;
        if(der[4]) der[4][i] = dpy * arg_y;
        return val*height;
    }
    inline double operator ()(const double &x, const double &y){
        double arg_x = (x-pos_x)/width_x;
        double arg2_x = arg_x*arg_x;
        double arg_y = (y-pos_y)/width_y;
        double arg2_y = arg_y*arg_y;
        double arg2 = arg2_x+arg2_y;
        double val = exp(-gauss_width_cof*arg2);
        return height * val;
    }
private:
    double height;
    double pos_x;
    double pos_y;
    double width_x;
    double width_y;
};
#endif
