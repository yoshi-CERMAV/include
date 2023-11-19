/*
 *  rotator.h
 *  ghemical
 *
 *  Created by Yoshiharu Nishiyama on 05/12/13.
 *  Copyright 2005 Cermav. All rights reserved.
 *
 */
 
#ifndef ROTATOR_H
#define ROTATOR_H 
#include <iostream>
#include <math.h>
#include <string.h>
using namespace std;
template <class T>
class rotator
{
protected:

   T *matrix;
   
public:
   rotator(){
     matrix = new T[9];
	 init();
   }
  
   void swap(){ //swapping x and z coordinate
      T t= matrix[2]; matrix[2] = matrix[0]; matrix[0] = t;
      t= matrix[5]; matrix[5] = matrix[3]; matrix[3] = t;
      t= matrix[8]; matrix[8] = matrix[6]; matrix[6] = t;
   } 
   void dump(){
   cout << "matrix\n"<<matrix[0]<<" "<<matrix[1]<<" "<<matrix[2]<<endl;
   cout << matrix[3]<<" "<<matrix[4]<<" "<<matrix[5]<<endl;
   cout << matrix[6]<<" "<<matrix[7]<<" "<<matrix[8]<<endl;
   }
   void get_view(){
       cout <<matrix[0] <<", "<<matrix[1] <<", "<<matrix[2]<<",\\"<<endl;
       cout <<matrix[6] <<", "<<matrix[7] <<", "<<matrix[8]<<",\\"<<endl;
       cout <<matrix[3] <<", "<<matrix[4] <<", "<<matrix[5]<<",\\"<<endl;

   }
   
   void init()
   {
  	  for(int i = 0; i < 9; i++){ 
	    if(!(i%4)) matrix[i] = 1; 
	    else matrix[i] = 0;
     }
   }
   
   inline void rot(const T&x, T *a0, T *a1){      
      T c = cos(x);
	  T s = sin(x);
	  for(int i = 0; i < 3; i++){
	    T t2 = c* *(a1+i) - s* *(a0+i);
	    *(a0+i) =  c* *(a0+i) + s* *(a1+i);
	    *(a1+i) = t2;
	  }
  }

   void rotx(T x){
         rot(x, matrix+3, matrix+6);
   }
   void roty(T x){
         rot(x, matrix+6,  matrix);
   }
   void rotz(T x){
         rot(x, matrix, matrix+3);
   }
   void apply(T *x){
     T t0 = x[0]*matrix[0]+x[1]*matrix[1]+x[2]*matrix[2];
     T t1 = x[0]*matrix[3]+x[1]*matrix[4]+x[2]*matrix[5];
	 x[2] = x[0]*matrix[6]+x[1]*matrix[7]+x[2]*matrix[8];
	 x[0] = t0;
	 x[1] = t1;
   }
    void bring_on_x(T *x){
        init();
        double a0 = atan2(x[1], x[0]);
        rotz(a0);
        T temp[3];
        memcpy(temp, x, sizeof(T)*3);
        apply(temp);
        double a1 = atan2(temp[2], temp[0]);
        roty(-a1);
        memcpy(temp, x, sizeof(T)*3);
        apply(temp);
        cout << temp[0] <<" "<< temp[1] << " "<< temp[2] <<endl;
    }
    void bring_on_xy (T *x, T *y){
        init();
        double a0 = atan2(x[1], x[0]);
        rotz(a0);
        T temp[3];
        memcpy(temp, x, sizeof(T)*3);
        apply(temp);
        double a1 = atan2(temp[2], temp[0]);
        roty(-a1);
        memcpy(temp, y, sizeof(T)*3);
        apply(temp);
        double a2 = atan2(temp[2], temp[1]);
        rotx(a2);
        memcpy(temp, y, sizeof(T)*3);
        apply(temp);
        cout << temp[0] <<" "<<temp[1] <<" "<< temp[2]<<endl;
        memcpy(temp, x, sizeof(T)*3);
        apply(temp);
        cout << temp[0] <<" "<<temp[1] <<" "<< temp[2]<<endl;

    }
};




template<class T>
T angle(const T &x, const T &y)  // returns between -pi:pi
{
   return atan2(y, x);
}

template<class T>
void
getangle(const T b[6], T ang[3])
{
  rotator<T> rot;
  if(fabs(b[1])>1e-6||fabs(b[2]) > 1e-6){
  ang[0] = angle(b[2], b[1]);
  rot.rotx(-ang[0]);
  }
  else ang[0] = 0;
  T temp[3];

  memcpy(temp, b,sizeof(T)*3);
  rot.apply(temp);
  if(fabs(temp[2])>1e-6||fabs(temp[0]) > 1e-6){  
  ang[1] = angle(temp[2], temp[0]);
  rot.roty(-ang[1]);
  }
  else ang[1] = 0;

  memcpy(temp,&(b[3]), sizeof(T)*3);
  rot.apply(temp);
  if(fabs(temp[1])>1e-6||fabs(temp[0]) > 1e-6){  
  ang[2] = angle(temp[1], temp[0]);
  }else ang[2] = 0;
}


template<class T>
void
get_angle_zy(const T b0[3],  T ang[2])
{
  rotator<T> rot;
  if(fabs(b0[0])>1e-6||fabs(b0[1]) > 1e-6){
  ang[0] = angle(b0[0], b0[1]);
  rot.rotz(ang[0]);
  }
  else ang[0] = 0;
  T temp[3];

  memcpy(temp, b0,sizeof(T)*3);
//  cout <<"temp "<<temp[0]<<" "<<temp[1]<<" "<<temp[2]<<endl;  
  rot.apply(temp);
//  cout <<"temp "<<temp[0]<<" "<<temp[1]<<" "<<temp[2]<<endl;  
  if(fabs(temp[2])>1e-6||fabs(temp[0]) > 1e-6){  
  ang[1] = angle(temp[2], temp[0]);
  rot.roty(ang[1]);
  }
  else ang[1] = 0;
  memcpy(temp, b0,sizeof(T)*3);
  rot.apply(temp);
//  cout <<"temp "<<temp[0]<<" "<<temp[1]<<" "<<temp[2]<<endl;  
}

template <class T>
void 
get_angle_xy(const T b0[3], T ang[2])
{
   rotator <T> rot;
   ang[0] = atan2(b0[1], b0[2]);
   rot.rotx(ang[0]);
   T temp[3];
   memcpy(temp ,b0 ,sizeof(T)*3);
   rot.apply(temp);
   ang[1] = atan2(b0[0], b0[2]);

}





#endif

