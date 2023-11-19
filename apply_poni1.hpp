//
//  apply_poni1.hpp
//  D2AM2023
//
//  Created by yoshi on 09/08/2023.
//

#ifndef apply_poni1_hpp
#define apply_poni1_hpp
#include <hdf5.h>
#include <hdf5_hl.h>
#include <H5Cpp.h>
#include <iostream>
#include <fstream>
#include <Accelerate/Accelerate.h>
#include <map>
#include <stdio.h>
#include <cmath>
#include "rotator.h"
const double four_pi = M_PI * 4;

using namespace H5;
using namespace std;

class Comparator
{
public:
    Comparator(float *d){dat = d;}
    bool operator() (size_t i, size_t j){return dat[i*2] < dat[j*2];}
private:
    float *dat;
};

class detector
{
public:
    detector(){}
    void init(const char H5[]){
        H5File file(H5, H5F_ACC_RDONLY);
        DataSet dataset = file.openDataSet("/entry_0000/pyFAI/Detector/pixel_corners");
        DataSpace dataspace = dataset.getSpace();
        
        int rank = dataspace.getSimpleExtentNdims() ;
        hsize_t dims_out[rank];
        const int ndims = dataspace.getSimpleExtentDims( dims_out, NULL);
        pixel_corners_dim = 1;
        cout << "rank is "<< rank <<" "<<ndims<<endl;
        for(int i = 0; i < rank; i++) {
            cout << dims_out[i]<<endl;
            pixel_corners_dim *= dims_out[i];
        }
        dim0 = dims_out[1];
        dim1 = dims_out[0];
        image_size = dims_out[0] * dims_out[1];
        pixel_corners = new float [pixel_corners_dim];
        center_pos = new float[image_size*3];
        cout << "pixel_corners_dim" <<pixel_corners_dim <<endl;
        cout << "image_size = " << image_size<<endl;
        cout << "dim0, dim1 = " << dim0<<", "<<dim1<<endl;
        dataset.read(pixel_corners, PredType::IEEE_F32LE, dataspace, dataspace);
        
        float *ptr1 = center_pos;
        float *ptr = pixel_corners;
        for(int i = 0; i < image_size; i++, ptr1+=3, ptr +=12){
            //       cout << i << endl;
            for(int j = 0; j < 3; j++){
                ptr1[j] = ptr[j] + ptr[j+3] + ptr[j+6] + ptr[j+9];
                ptr1[j] *= 0.25;
                //            cout << ptr1[j]<<endl;
            }
        }
        get_min_max(min_x, max_x ,min_y, max_y);
        cout << "min_x = "<<min_x <<endl;
        cout << "max_x = "<<max_x <<endl;
        cout << "min_y = "<<min_y <<endl;
        cout << "max_y = "<<max_y <<endl;
        float shape = (max_y-min_y)/(max_x-min_x);
        cout << "shape = "<<shape<<endl;
    }
    detector(const char H5[])
    {
        init(H5);
    }
    void get_min_max(float &min_x, float &max_x, float &min_y, float &max_y)
    {
        min_x = 1e5;
        max_x = -1e5;
        min_y = 1e5;
        max_y = -1e5;
        float *ptr = center_pos;
        for(int i = 0; i < image_size; i++, ptr+=3){
            if(ptr[2] < min_x ) min_x = ptr[2];
            if(ptr[2] > max_x ) max_x = ptr[2];
            if(ptr[1] < min_y ) min_y = ptr[1];
            if(ptr[1] > max_y ) max_y = ptr[1];
        }
    }
    
    int find_pos(float x, float y){
        float *ptr = center_pos;
        float min = 1e5;
        int pos;
     //   cout << x<<" "<<y << endl;

        for(int i = 0; i < image_size; i++, ptr+=3){
            float dx, dy;
            if((dx=(x-ptr[2])) > step) continue;
            if((dy=(y-ptr[1])) > step) continue;
            float dist = dx*dx + dy*dy;
            if (dist < min){
                min = dist;
                pos = (ptr-center_pos)/2;
      //          cout << x<<" "<<y << " "<< dist <<" "<<pos<<endl;

            }
        }
        if (min < step) return pos;
        else return -1;
    }
    void init_quick_view(int width, int height){
    //    height = int ( width * shape);
        step = ( max_x - min_x) / width;
        cout << " step = "<<step<<endl;
        view_size = width * height;
        pixel_pos = new int[view_size];
        cout <<" allocated pixel_pos"<<endl;
        int count = 0;
        int *iptr = pixel_pos;
        float *dist = new float[view_size];
        for(int i = 0; i < view_size; i++) pixel_pos[i] = -1;
        cout <<" allocated dist"<<endl;
        for(int i = 0; i < view_size; i++) dist[i] = 1e5;
        cout <<"initialize dist"<<endl;
        for(int i = 0; i < image_size; i++){
//            if (!(i%1000)) cout << i << endl;
            float x = center_pos[i*3+2];
            float y = center_pos[i*3+1];
            float x1 = x/step;
            float y1 = y/step;
            int xi = int(x1+0.5);
            int yi = int(y1+0.5);
            float xs = x1-xi;
            float ys = y1-yi;
            if (xi < 0 ) continue; if (xi > (width-2)) continue;
            if (yi < 0 ) continue; if (yi > (height-2)) continue;
            int pos = xi + width * yi;
            float a = xs*xs + ys*ys ;
            if(a < dist[pos]) {
                dist[pos] = a;
                iptr[pos] = i;
            }
        }
        cout <<"initialized"<<endl;
    }
    
    void make_image1(unsigned char *image, int *data, float x0, float x1){
        double scale = 255./(x1-x0);
        make_image(image, data, x0, scale);
    }
    void make_image(unsigned char *image, int *data, float x0, float scale){
     //   cout <<"view size "<<view_size<<endl;
        if(!image)cerr<<"make_image error: no image space";
        if(!data)cerr<<"make_image error: no data space";
        for(int i = 0; i < view_size; i++){
            int ipos = pixel_pos[i];
            if(ipos<0) {
                image[i] = 0;}
            else{
                float a = (data[ipos] - x0) * scale;
                if (a<0) a= 0;
                if (a>255)a=255;
                image[i] = (char)a;
            }
        }
    }
    void dump_pixel_corners(){
        float *ptr = pixel_corners;
        for(int i = 0; i < 4; i++, ptr+=3){
            cout << ptr[0] << " "<<ptr[1]<<" "<<ptr[2]<<endl;
        }
        cout <<endl;
        ptr = pixel_corners;
        for(int i = 0; i < 5; i++, ptr+=12){
            cout << ptr[0] << " "<<ptr[1]<<" "<<ptr[2]<<endl;
        }
        cout << endl;
        ptr = pixel_corners+12*dim0*120;
        for(int i = 0; i < 5; i++, ptr+=12){
            cout << ptr[0] << " "<<ptr[1]<<" "<<ptr[2]<<endl;
        }

    }
protected:
    size_t dim0, dim1;
    size_t image_size;
    size_t pixel_corners_dim;
    float *pixel_corners = NULL;
    float *center_pos = NULL;
    int *pixel_pos;
    float min_x, max_x, min_y, max_y;
    float shape;
    float step;
    size_t view_size;
    
};

class apply_poni:public detector
{
public:
    apply_poni(const char H5[], const char poni_filename[]):detector(H5){
        allocate();
        read_poni(poni_filename);
        dump_pixel_corners();
        init_map();
    }
    apply_poni(const char filename[]):detector(){
        ifstream fi(filename);
        char h5name[255];
        char poni_name[255];
        fi.getline(h5name, 255);
        cout << h5name << endl;
        detector::init(h5name);
        fi.getline(poni_name,255);
        cout << poni_name <<endl;
        allocate();
        read_poni(poni_name);
        dump_pixel_corners();
        init_map();

    }
    
    void read_poni(const char poni_filename[])
    {
        char key[255];
        double value;
        char *buffer = new char [1024];
        
        // reading ponifiel
        ifstream ponifile(poni_filename);
        
        while(ponifile.getline(buffer, 1024)){
            sscanf(buffer,"%[^:]:%lf",key, &value);
            poni[key]=value;
            cout << key <<": "<< value <<endl;
        }
        
        cout <<": "<<  poni["Poni1"]<<endl;
        r.init();
        r.roty(-poni["Rot1"]);
        r.rotx(-poni["Rot2"]);
        r.swap(); // swapping x and z;
        center[0] = poni["Distance"];
        center[1] = -poni["Poni1"];
        center[2] = -poni["Poni2"];
        cout << "center "<<center[0] <<" "<<center[1]<<" "<<center[2]<<endl;
        wavelength = poni["Wavelength"];
    }
    void allocate(){
        qb = new float[8*image_size];
        pixel_index = new size_t[image_size];
        center_qb = new float[image_size * 2];
        center_corr = new float[image_size];
        mask = new char[image_size];
        for(int i = 0; i != image_size; i++) mask[i] = 0;
    }
    
    void set_weight(double *flat){
        if(!weight) weight = new float[image_size];
        for(int i = 0; i < image_size; i++){
            if (flat[i]) weight[i] = 1./flat[i];
            if(flat[i] < 0.5 || flat[i] > 2) mask[i] = 1;
        }
    }
    
    double get_q(int x, int y){
        size_t pos = y*dim0 + x;
        return center_qb[pos*2];
    }
    
    void init_map()
    {
        float *ptr = pixel_corners;
        float *ptr_qb= qb;
        scale = four_pi / wavelength *1e-10;
        int num_point = image_size * 4;
        for(int i = 0; i < num_point ; i++, ptr+=3, ptr_qb+=2){
            add(ptr, center);
            r.apply(ptr);
            double r2 = sqrt(ptr[0] * ptr[0] + ptr[1] * ptr[1]);
            double theta = 0.5 * atan2(r2, ptr[2]);
            ptr_qb[0] = scale * sin(theta);
            ptr_qb[1] = atan2(ptr[0], ptr[1]);
        }
        cout << "checking if on the same side"<<endl;
        float M_PI2 = 2*M_PI;
        for(int i = 0; i < image_size; i++){
            float *ptr = qb+i*8+1;
            for(int j = 2; j < 8; j+=2){
                if(fabs(ptr[0]-ptr[j])>M_PI){
                    ptr[j] = (ptr[j] >0 )?ptr[j]-M_PI2 :ptr[j]+M_PI2;
                }
            }
        }
        cout << "done"<<endl;
        // calculate center of pixels
        cout << "calculating correction"<<endl;
        //calculate correction factors
        ptr = center_pos;
        ptr_qb = center_qb;
        float norm[3];
        norm[0] = 1; norm[1] = 0; norm[2] = 0;
        r.apply(norm);
        for(int i = 0; i < image_size; i++, ptr+=3, ptr_qb+=2){
            double r2 = sqrt(ptr[0] * ptr[0] + ptr[1] * ptr[1]);
            double Theta = atan2(r2, ptr[2]);
            ptr_qb[0] = scale * sin(0.5*Theta);
            double beta =atan2(ptr[0], ptr[1]);
            ptr_qb[1] =beta;
            double cT = cos(Theta);
            double sT = sin(Theta);
            cT*=cT; sT*=sT;
            float corr = 0.5*(1 + cT + 1*cos(beta * 2) * sT); // phase 90˚ diff from -cos..
            corr = 1./corr;
            double r3 =ptr[0] * ptr[0] + ptr[1] * ptr[1] + ptr[2]*ptr[2];
            corr *= fabs(r3*sqrt(r3)/scalar(ptr, norm));
            center_corr[i] = corr;
        }
        cout <<"going to sort "<<endl;
        // sort as function of q.
        Comparator comp(center_qb);
        for(int i = 0; i!= image_size; i++) pixel_index[i]=i;
        sort(pixel_index, pixel_index+image_size, comp);
        cout <<"sorted"<<endl;
        auto [q, r] = std::div((int)pixel_index[0], 578);
        cout << "q r "<< q<<" "<<r<<" "<<center_qb[pixel_index[0]*2]<<endl;
        auto [q1, r1] = std::div((int)pixel_index[image_size-1], 578);
        cout << "q r "<< q1<<" "<<r1<<" "<<center_qb[pixel_index[image_size-1]*2]<<endl;
    }
    void add(float *a, float *b)
    {
        a[0] += b[0];
        a[1] += b[1];
        a[2] += b[2];
    }
    float sq(float *x)
    {
        return x[0] * x[0] + x[1]*x[1] + x[2]*x[2];
    }
    float scalar(float *a, float *b)
    {
        return a[0]*b[0] + a[1]*b[1] + a[2]*b[2];
    }
    void write(ofstream &fo, double data){
        fo.write(reinterpret_cast<char *>(&data), sizeof(double));
    }
    void write(ofstream &fo, int data){
        fo.write(reinterpret_cast<char *>(&data), sizeof(int));
    }

    void extract_intensity(double qmin, double qmax, int *dat, double *flat, char filename[])
    {
        ofstream fo(filename);
        int i=0, pix;
        double q;
        while(center_qb[pixel_index[i]*2] < qmin) i++;
        do{
            pix = pixel_index[i];
            q =center_qb[pix*2];
            write(fo, q);
            write(fo, center_qb[pix*2+1]);
            write(fo, dat[pix]*flat[pix]);
            i++;
        }while(q <qmax);
    }
    
    void extract_intensity1(double qmin, double qmax, int *dat, double *flat, char filename[])
    {
        ofstream fo(filename);
        int i=0, pix;
        double q;
        while(center_qb[pixel_index[i]*2] < qmin) i++;
        do{
            pix = pixel_index[i];
            q =center_qb[pix*2];
            write(fo, q);
            write(fo, center_qb[pix*2+1]);
            write(fo, dat[pix]*flat[pix]*center_corr[pix]);
            i++;
        }while(q <qmax);
    }
    void add_mask(const char *mask1){
        for(int i = 0; i != image_size; i++) if(mask1[i]) mask[i] = 1;
    }
    int mask_line(int n){
        char *ptr = mask + n*dim0;
        for(int i = 0; i != dim0; i++ ) ptr[i] = 1;
        return 0;
    }
    void integrate(double qmin, double qmax, double qstep,
                   int *bg, int *dat, double scale, char filename[])
    {
        ofstream fo(filename);
        double span = qstep *0.5;
        int kmax = (qmax-qmin)/qstep+1;
        double q0 = qmin - span;
        double q1 = qmin + span;
        int i0 = 0;
        for(int k = 0; k < kmax; k++){
            while (center_qb[pixel_index[i0]*2] < q0 && i0< image_size) i0++;
            int i1 = i0;
            while (center_qb[pixel_index[i1]*2] < q1 && i1 < image_size) i1++;
//            cout << i1 <<" "<<i0<<endl;
            size_t len = i1-i0;
            size_t *pix = pixel_index + i0;
            double sum = 0;
            double sum2 = 0;
            double count = 0;
            for(int i = 0; i!=len; i++){
                size_t p = pix[i];
                if(mask[p])continue;
                double w = weight[p];
                double my_bg =bg[p];
                double my_data = dat[p];
                double v = my_data - scale * my_bg;
                v*= center_corr[p];
                sum+=v;
                sum2+=v*v/w;
                count+=w;
            }
            if(count){
                sum/=count;
                sum2/=count;
                double stdev = sqrt(sum2-sum*sum);
//                cout<< q0 + span<<" "<<sum/scale <<" "<<stdev/scale <<" "<<count <<endl;

                write(fo, q0+span);
                write(fo, sum);
                write(fo, stdev);
                write(fo, count);
            }
            q0 += qstep;
            q1 += qstep;
        }
    }
    
    void integrate(double qmin, double qmax, double qstep,
                   int *bg, int *dat, double *flat, double scale, char filename[])
    {
        double inv_scale = 1./scale;
        ofstream fo(filename);
        double span = qstep *0.5;
        int kmax = (qmax-qmin)/qstep+1;
        double q0 = qmin - span;
        double q1 = qmin + span;
        int i0 = 0;
        for(int k = 0; k < kmax; k++){
            while (center_qb[pixel_index[i0]*2] < q0 && i0< image_size) i0++;
            int i1 = i0;
            while (center_qb[pixel_index[i1]*2] < q1 && i1 < image_size) i1++;
            size_t len = i1-i0;
            size_t *pix = pixel_index + i0;
            double sum = 0;
            double sum2 = 0;
            int count = 0;
            for(int i = 0; i!=len; i++){
                size_t p = pix[i];
                double f = flat[p];
                if(f< 0.5 || f > 2 )continue;
                double my_bg =bg[p];
                double my_data = dat[p];
                double v = my_data * inv_scale - my_bg;
                v*= f;
                v*= center_corr[p];
                sum+=v;
                sum2+=v*v;
                count++;
            }
            if(count){
                sum/=count;
                sum2/=count;
                double stdev = sqrt(sum2-sum*sum);
                
                write(fo, q0+span);
                write(fo, sum);
                write(fo, stdev);
                write(fo, count);
                //                fo << q0 + span<<" "<<sum/scale <<" "<<stdev/scale <<" "<<count <<endl;
            }
            q0 += qstep;
            q1 += qstep;
        }
    }
    // fit against "int" data
    double fit1(double q0, double q1, int *dat, int *bg, double *flat, double *flat_err){
        int i0 = 0;
        while (center_qb[pixel_index[i0]*2] < q0) i0++;
        int i1 = i0;
        while (center_qb[pixel_index[i1]*2] < q1) i1++;
        int len = i1-i0;
        double *A = new double [len * 4];
        int *pivot = new int[5];
        double *y = A + len *3;
        size_t *pix = pixel_index + i0;
        double *A1 = A+len;
        double *A2 = A1+len;
        cout <<" len "<<len<<endl;;
        ofstream fo("/Users/yoshi/temp/dump.txt");
        for(int i = 0; i!=len; i++){
            size_t p = pix[i];
            double flat_Value = flat[p];
            double flat_Err = flat_err[p];

            if(isnan(flat_Value) || isnan(flat_Err)||flat_Value < .90 || flat_Value > 1.1  || mask[p]){
                A[i] = y[i] = A1[i] = A2[i] = 0;
       //         cout << "isnan "<<i<<endl;
                continue;
            }else{
                double my_bg =bg[p];
                double my_data = dat[p];
//                double my_sum = my_bg + my_data;
                
                A[i] = my_bg;
                y[i] = my_data;
                fo << i <<" "<<center_qb[p*2] <<" "<<center_qb[p*2+1]<<" "<<A[i]<<" "<<y[i];
                A1[i] = 1./flat_Value;
                A2[i] = center_qb[p*2]/flat_Value;
#define AAA
#ifdef AAA
                double sigma = 1./sqrt(flat_Value);//(my_sum) * flat_Err + sqrt(my_sum) * flat_Value;
                if (sigma < 0){
                    A[i] = 0;
                    y[i] = 0;
                    A1[i] = 0;
                    A2[i] = 0;
                }else{
 //                   sigma = 1./sigma;
                    A[i] *= sigma;
                    y[i] *= sigma;
                    A1[i] *= sigma;
                    A2[i] *= sigma;
                    //              }
                }
                fo <<" "<<sigma<<endl;
#endif
            }
        }
        cout <<"dgelsy "<<endl;
        int m = 2;
        int one = 1;
        double rcond = 1e-32;
        int rank1;
        int info;
        double  iwork=-1;
        int lwork = -1;
        dgelsy_(&len, &m, &one, A, &len,
                y, &len, pivot, &rcond,
                &rank1, &iwork, &lwork, &info);
        cout << iwork <<endl;
        lwork = iwork;
        double *work = new double [lwork];
        dgelsy_(&len, &m, &one, A, &len,
                y, &len, pivot, &rcond,
                &rank1, work, &lwork, &info);
        cout << y[0] <<" "<< y[1]<<" "<<y[2]<<" "<< rank1<<endl;
        return y[0];
    }
    
    // fit against cleaned "double" data.
    double fit(double q0, double q1, double *dat, double *bg, double *flat, double *flat_err){
        int i0 = 0;
        while (center_qb[pixel_index[i0]*2] < q0) i0++;
        int i1 = i0;
        while (center_qb[pixel_index[i1]*2] < q1) i1++;
        int len = i1-i0;
        double *A = new double [len * 4];
        int *pivot = new int[5];
        double *y = A + len *3;
        size_t *pix = pixel_index + i0;
        double *A1 = A+len;
        double *A2 = A1+len;
        cout <<" len "<<len<<endl;;
        ofstream fo("/Users/yoshi/temp/dump.txt");
        for(int i = 0; i!=len; i++){
            size_t p = pix[i];
            double flat_Value = flat[p];
            double flat_Err = flat_err[p];

            if(isnan(flat_Value) || isnan(flat_Err)||flat_Value < .90 || flat_Value > 1.1){
                A[i] = y[i] = A1[i] = A2[i] = 0;
       //         cout << "isnan "<<i<<endl;
                continue;
            }else{
                double my_bg =bg[p];
                double my_data = dat[p];
                double my_sum = my_bg + my_data;
                if(isnan(my_bg) || isnan(my_data)){
                    A[i] = y[i] = A1[i] = A2[i] = 0;
                    continue;
                }
                A[i] = my_bg * flat_Value;
                y[i] = my_data * flat_Value;
                fo << i <<" "<<center_qb[p*2] <<" "<<center_qb[p*2+1]<<" "<<A[i]<<" "<<y[i];
                A1[i] = 1;
                A2[i] = center_qb[p*2];
//#define AAA
#ifdef AAA
                double sigma = (my_sum) * flat_Err + sqrt(my_sum) * flat_Value;
                if (sigma < 0){
                    A[i] = 0;
                    y[i] = 0;
                    A1[i] = 0;
                    A2[i] = 0;
                }else{
                    sigma = 1./sigma;
                    A[i] *= sigma;
                    y[i] *= sigma;
                    A1[i] *= sigma;
                    A2[i] *= sigma;
                    //              }
                }
                fo <<" "<<sigma<<endl;
#endif
            }
        }
        cout <<"dgelsy "<<endl;
        int m = 2;
        int one = 1;
        double rcond = 1e-32;
        int rank1;
        int info;
        double  iwork=-1;
        int lwork = -1;
        dgelsy_(&len, &m, &one, A, &len,
                y, &len, pivot, &rcond,
                &rank1, &iwork, &lwork, &info);
        cout << iwork <<endl;
        lwork = iwork;
        double *work = new double [lwork];
        dgelsy_(&len, &m, &one, A, &len,
                y, &len, pivot, &rcond,
                &rank1, work, &lwork, &info);
        cout << y[0] <<" "<< y[1]<<" "<<y[2]<<" "<< rank1<<endl;
        return y[0];
    }
protected:
    float *qb = NULL;
    float *center_qb = NULL;
    float *center_corr = NULL;
    float *weight=NULL;
    size_t *pixel_index = NULL;
    char *mask;
    map<string, double> poni;
    rotator<float> r;
    float center[3];
    double wavelength;
    float scale;
    
};
#endif /* apply_poni1_hpp */
