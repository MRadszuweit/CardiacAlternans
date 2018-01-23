
#include "fft.h"

    
complex* FFT_data = NULL;
complex* FFT_trans = NULL;    
int Data_size;
int FFT_mode;
double FFT_sign;

complex cmp_sum(complex a, complex b){
	complex res;
	res.x = a.x+b.x;
	res.y = a.y+b.y;
	return res;
}

complex cmp_diff(complex a, complex b){
	complex res;
	res.x = a.x-b.x;
	res.y = a.y-b.y;
	return res;
}

complex cmp_prod(complex a, complex b){
	complex res;
	res.x = a.x*b.x-a.y*b.y;
	res.y = a.x*b.y+a.y*b.x;
	return res;
}

double cmp_sqrabs(complex z){
	return z.x*z.x+z.y*z.y;
}

void phasemult(complex* z,double phase){
        double a = z->x*cos(phase)-z->y*sin(phase);
        double b = z->x*sin(phase)+z->y*cos(phase);
        z->x = a;
        z->y = b;
}

complex* cmp_expand(double* data,int size){
	int i;
	complex* Res = (complex*)malloc(size*sizeof(complex));
	for (i=0;i<size;i++){
		Res[i].x = data[i];
		Res[i].y = 0;
	}
	return Res;
}

complex* cmp_cpy(complex* Z,int size){
	int i;
	complex* Res = (complex*)malloc(size*sizeof(complex));
	for (i=0;i<size;i++){
		Res[i].x = Z[i].x;
		Res[i].y = Z[i].y;
	}
	return Res;
}
    
void init_FastFourier1D(complex* data,int size){
        FFT_data = data;
        Data_size = size;
        FFT_mode = ANALYSE;
        FFT_sign = -1;
    }
    
void set_FFT_mode(int mode){
        FFT_mode = mode;
        if (mode==ANALYSE) FFT_sign = -1; else FFT_sign = 1;
}
    
index2D map(index2D c,int even){
        index2D res;
        if (even){
			res.i = c.i;
			res.j = 2*c.j;
		}
        else{
			res.i = c.i+c.j;
			res.j = 2*c.j;
		}
        return res;
}

complex* iterate(index2D c,int len){
        complex* res = (complex*)malloc(len*sizeof(complex));
        if (len==2){
            complex p1 = FFT_data[c.i];
            complex p2 = FFT_data[c.i+c.j];
            res[0] = cmp_sum(p1,p2);
            res[1] = cmp_diff(p1,p2);           
        }        
        else{
			int i;
            int m = len/2;
            complex* f0 = iterate(map(c,EVEN),m);
            complex* f1 = iterate(map(c,ODD),m);
            for (i=0;i<m;i++){
				phasemult(&(f1[i]),(double)2.*FFT_sign*M_PI*i/len);
				res[i] = cmp_sum(f0[i],f1[i]);         
				res[i+m] = cmp_diff(f0[i],f1[i]);
            }
            free(f0);
            free(f1);
        }
        return res;
}
    
void normalize(){
        if (FFT_trans!=NULL){
			int i;
            for (i=0;i<Data_size;i++){
				FFT_trans[i].x /= (double)Data_size;
				FFT_trans[i].y /= (double)Data_size;
			}
        }
}
    
void fastfourier(){
        if (FFT_data!=NULL){
			index2D c;
			c.i = 0;
			c.j = 1;
            //printf("Fast-Fourier-Transformation 1D\n");
            if (FFT_trans!=NULL) free(FFT_trans);
            FFT_trans = iterate(c,Data_size);
            if (FFT_mode==ANALYSE){normalize();}
        }
}

complex* get_FFT_trans(){
	return FFT_trans;
}

double* getstrfunction(){
        if (FFT_trans!=NULL){
			int i;
            double* Res = (double*)malloc(Data_size*sizeof(double));           
            for (i=0;i<Data_size;i++) Res[i] = cmp_sqrabs(FFT_trans[i]);
            return Res;
        }
        else return NULL;
}
