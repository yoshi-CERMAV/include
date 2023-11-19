//
//  mib-read.h
//  
//
//  Created by Yoshiharu Nishiyama on 15/11/2023.
//

#ifndef mib_read_h
#define mib_read_h
#include <iostream>
#include <fstream>
#include <sys/stat.h>
#include <Accelerate/Accelerate.h>
#include <Savitzky_Golay_2d.h>

using namespace std;
template<class T>
void center_of_gravity(T *data, int w, int h, int line_len, double &x, double &y)
{
   double sum= 0;
   double sumx = 0;
   double sumy = 0;
   T *ptr = data;
   for(int j = 0; j < h; j++, ptr += line_len){
      for(int i = 0; i < w; i++){
          T value =ptr[i];
          sum += value;
          sumx += value * i;
          sumy += value * j;
      }
   }
   x = sumx/sum;
   y = sumy/sum;
}

class mib
{
    const int line_len = 514;
    const int num_row = 514;
    const int image_size = 514*514;
    const int head_size = 768;
    const int frame_size = image_size * 2;
    const int frame_step = frame_size + head_size;
    size_t max_frame;
    ifstream fi;
    unsigned short *data;
    double *bg;
    double *d_data, *o_data;
    double *A, *B, *work;
    double scale;
    double max_value;
    double cx, cy;

    void swap()
    {
       char *cptr = reinterpret_cast<char *>(data);
       for(int i = 0; i < image_size ;i++, cptr+=2){
         char temp = *cptr;
         cptr[0] = cptr[1];
         cptr[1] = temp;
       }
    }
    void scale_hline(const int i)
    {
        unsigned short *ptr = data+i*line_len;
        for(int j = 0; j < line_len; j++){
            *(ptr++) *= scale;
        }
    }
    void scale_vline(const int i)
     {
        unsigned short *ptr = data+i;
        for(int j = 0; j < num_row; j++){
            *(ptr+= line_len) *= scale;
         }
     }
    void  copy_hline(int out_i, int in_i)
    {
        unsigned short *ptr = data + in_i * line_len;
        unsigned short *ptr1 = data + out_i * line_len;
        for(int i = 0; i < line_len; i++) *(ptr1++) = *(ptr++);
    }
    void copy_vline(int out_i, int in_i)
     {
        unsigned short *ptr = data + in_i;
        unsigned short *ptr1 = data + out_i;
        for(int i = 0; i < num_row; i++){
            *(ptr1+=line_len) = *(ptr+=line_len);
        }
     }
public:
    mib(char filename[], double in_scale){
        struct stat st;
        if(stat(filename, &st) != 0) {
            cerr << "file does not exist";
        }
        fi.open(filename);
        max_frame = st.st_size / frame_step;
        cout << max_frame<<endl;
        data = new unsigned short[image_size];
        d_data = new double[image_size];
        o_data = new double[image_size];

        A = new double [image_size];
        B = new double [image_size];
        work = new double[image_size];
        scale = in_scale;
    }
    
    void read(int frame){
        if (frame >= max_frame) return;
        fi.seekg(frame_step*frame + head_size, ios_base::beg);
        fi.read(reinterpret_cast<char *>(data), frame_size);
        swap();
        scale_hline(255);
        scale_hline(258);
        copy_hline(256, 255);
        copy_hline(257, 258);
        scale_vline(255);
        scale_vline(258);
        copy_vline(256, 255);
        copy_vline(257, 258);
    }
    int copy(unsigned char *cptr)
    {
        unsigned short *sptr = data;
        float ratio = 255./max_value;
        float temp;
         for(int i = 0; i < image_size; i++, sptr++){
            *(cptr++) = (temp = *(sptr) *ratio) > 255? 255: temp ;
         }
        return 0;
    }
    
    int copy(unsigned char *cptr, double *dat)
    {
        double *ptr = dat;
        double ratio = 255./max_value;
        double temp;
        for(int i = 0; i < image_size; i++, ptr++){
            *(cptr++) = (temp = *(ptr) *ratio) > 255? 255: temp ;
        }
        return 0;
    }
    int copy_o(unsigned char *cptr)
    {
        return copy(cptr, o_data);
    }
    int copy_d(unsigned char *cptr)
    {
        return copy(cptr, d_data);
    }
    
    void find_center(double &x, double &y, int size){
        int ix = x-size;
        int iy = y-size;
        if(ix <0){cerr<< "size too big- x"<<endl; return;}
        if(iy <0){cerr<< "size too big- y"<<endl;return;}
        int box = size*2+1;
        if(ix+box >= line_len){cerr<< "size too big+ x"<<endl; return;}
        if(iy+box >= num_row){cerr<< "size too big+ y"<<endl; return;}
        center_of_gravity(data+iy*line_len+ix, box, box, line_len, x, y);
        x+=ix;
        y+=iy;
    }
    
    void find_center(){
        double x, y;
        center_of_gravity(data, 514, 514, 514, x, y);
        find_center(x, y, 50);
        find_center(x, y, 20);
        find_center(x, y, 10);
        cx = x;
        cy = y;
    }
    
    void bg_subtract(double C){
        double *ptr = d_data;
        unsigned short *sptr = data;
        for(int j = 0; j < num_row; j++){
            for(int i = 0; i < line_len; i++, ptr++, sptr++){
                double x=  i-cx;
                double y = j-cy;
                double r2 = x*x + y*y;
                if(*sptr < 4075)
                    *ptr = (*sptr) - C/r2;
                else *ptr = -1;
            }
        }
    }
    double bg_scale(double radius)
    {
        int i0 = cx - radius;
        int i1 = cx + radius+1;
        int j0 = cy - radius;
        int j1 = cy + radius+1;
        double radius2 = radius*radius;
        double *ptrA = A;
        double *ptrB = B;
        int count = 0;
        for(int j = j0; j < j1; j++){
            unsigned short *ptr = data + j*line_len;
            for(int i = i0; i < i1; i++){
                double x=  i-cx;
                double y = j-cy;
                double r2 = x*x + y*y;
                if(ptr[i]< 4000 && r2 < radius2 ){
                    *(ptrA++) = 1./r2;
                    *(ptrB++) = ptr[i];
                    count++;
                }
            }
        }
        char trans = 'N';
        int one = 1;
        int info;
        int lwork = image_size;
        cout <<"count "<<count << endl;
        if(count){
            dgels_(&trans, &count, &one, &one, A, &count, B, &count, work, &lwork, &info);
            cout << cx <<" "<<cy <<" "<<B[0]<<endl;
        }
        return B[0];
    }
    void apply(Filter2d &filter){filter.apply(d_data, o_data, 514 * 514); cout << "applied "<<endl;
  //      ofstream fo("d_data");
  //      fo.write(reinterpret_cast<char *> (d_data), image_size *sizeof(double) );
  //      ofstream fo1("o_data");
  //      fo1.write(reinterpret_cast<char *> (o_data), image_size *sizeof(double) );

    }
    void set_center(double in_cx, double in_cy){cx = in_cx; cy = in_cy;}
    void set_max(double in_max){max_value = in_max;}
    void set_scale(double in_scale){scale = in_scale;}
    unsigned short *data_ptr(){return data;}
    double *filtered_data(){return o_data;}
    double *double_data(){return d_data;}
};

#endif /* mib_read_h */
