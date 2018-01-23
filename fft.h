
#ifndef FFT_H
#define FFT_H
  
#include <stdio.h>
#include <math.h>
#include <stdlib.h>


typedef struct COMPLEX{
	double x;
	double y;
}complex;
 
typedef struct INDEX2D{
	int i;
	int j;
}index2D; 
 
#define ANALYSE 1
#define SYNTHESE 0 
#define EVEN 1
#define ODD 0

complex* cmp_expand(double* data,int size);
complex* cmp_cpy(complex* Z,int size);
void init_FastFourier1D(complex* data,int size);
void set_FFT_mode(int mode);
void fastfourier();
complex* get_FFT_trans();
double* getstrfunction();

 #endif
 
