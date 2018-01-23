#include <stdio.h>
#include <math.h>
#include <malloc.h>
#include <sys/stat.h>
#include <string.h>
#include <dirent.h>
#include "../FEM/LinA/linear_algebra.h"
#include "fft.h"

//////// defines ///////////////////////////////////////////////

#define DYNAMIC 0
#define S1S2 1
#define S1S2_FINISH 2
#define FIXED 0
#define EQUILIBRIUM 1
#define INDEPENDENT 0
#define CONTINUATION 1


//////// global vars ///////////////////////////////////////////

int N;  // time steps
int M;  //space steps
int* T_size = NULL;  // # of Ts
int* DI_size = NULL;  // # of DIs
int* APD_size = NULL;  // # of APDs
int* CV_size = NULL;  // # of CVs
int* S2_index_DI = NULL;  // for S1S2 identification
int* S2_index_APD = NULL;
int* S2_index_CV = NULL;
int* S2_index_T = NULL;
int pacing_protocol;  // dynamic or S1S2
int pacing_width;  // x-width of the pacing pulse
 // bounds for not taking CV
int equilibrium_flag;  // set 1 if equilibrium reached
int block_flag;		// set 1 if block encountered
int glob_flag;		// set 1 if global excitation
int S1S2_counter;
int finish_mode;
int S1S2_index_APD;
int S1S2_index_DI;

double pacing_Q;  // charge of pacing pulse
double t_pace;  // pacing countdown

double P_Ta;
double P_a;

double L; // spatial length
double total_time; // sim time
double x;  // space
double t;  // time
double dx;  // space step
double dt;  // time step


double u_T;
double tau_T;
double tau_T2;

double u_o;
double u_u;
double theta_v;
double theta_w;
double theta_v_m;
double theta_o;
double tau_v_1m;
double tau_v_2m;
double tau_v_p;
double tau_w_1m;
double tau_w_2m;
double tau_w_p;
double k_w;
double u_w;
double tau_fi;
double tau_o1;
double tau_o2;
double tau_so1;
double tau_so2;
double k_so;
double u_so;
double tau_s1;
double tau_s2;
double k_s;
double u_s;
double tau_si;
double tau_w_inf;
double w_inf;

double alpha;
double c0;
double a_c;
double tau_c;
double c50;
double rho0;
double rho1;
double bmax;
double Gamma;

double D;  // diffusion

double U_s;  // sa threshold constant
double G_s;  // stretch activation coupling
double E;  // elastic constant
double F_T;  // active tension
double C1;  // Mooney-Rivlin 1. invariant
double C2;  // Mooney-Rivlin 2. invariant
double a;  // correction factor for non-Dirichlet conditions
double E_o; // outer elastic modulus
double l_o; // outer width

double T_pace1;  // pacing period
double T_pace2;  // S1 pacing if S1S2-protocol

double eps_APD;  // (APD^(n+1)-APD^n)/APD^n equilibrium detection
double eps_Newton;  // epsilon for all Newton methods used
double thr_u;  // threshold for APD detection

double* u = NULL;  // normalized voltage
double* v = NULL;  // gating v
double* w = NULL;  // gating w
double* s = NULL;  // gating s
double* c = NULL;  // free calcium c
double* b = NULL;  // bound calcium b
double* z = NULL;  // tropomyosin z
double* Ta = NULL;  // active tension
double* F = NULL;  // deformation gradient
double* Old_u = NULL;  // voltage of t-dt

double* Last_crossing = NULL;  // for APD/DI measurements

double** T = NULL;  // period
double** DI = NULL;  // diastolic interval
double** APD = NULL;  // action potential duration
double** CV = NULL;  // conduction velocity

//////// save options //////////////////////////////////////////

FILE** Time_series = NULL;  // storage for time series
FILE* Spacetime_u = NULL;  // storage spacetime plot
char output_name[512];
char output_dir[512];
int* st_times;  // index in space
int* st_indices;  // index in time
int nb_times;  // number of time series to store
int nb_space;  // number of spatial plots to store
int data_interval;  // number of timesteps between each data storage
int space_interval;  // number of space step beetween each stored
int CV_bound;

//////// programm //////////////////////////////////////////////

// initialization //

double fu_init(double x){
	return 0;
}

double fv_init(double x){
	return 1.;
}

double fw_init(double x){
	return 1.;
}

double fs_init(double x){
	return 0;
}

double fc_init(double x){
	return c0;
}

double fb_init(double x){
	return c0/(c0+rho1*(1.-1./Gamma)/rho0);
}

double fz_init(double x){
	return 0;
}

double fT_init(double x){
	return 0;
}

double fF_init(double x){
	return 1.;
}

void set_space_points(int m){
	int i;
	//st_indices[0] = m/2;
	//st_indices[1] = (int)round((double)0.95*m);
	for (i=0;i<nb_space;i++){
		st_indices[i] = (int)round((double)m*(i+1.)/(nb_space+1.));
	}
}

void set_time_points(int n){
	int i;
	for (i=0;i<nb_times-1;i++){
		st_times[i] = (int)round((double)13045./dt+i*1000.);
		//(int)round((double)4.*n/5.+i*1000.);
		//printf("%d:%d\n",i,st_times[i]);
	}
	st_times[nb_times-1] = N-1;
}

void init_param(){
	
	data_interval = 100;
	space_interval = 10; //2;
	
	
	//////// discretization
	dt = 0.01;  // 0.1 standart, 0.01 hier wichtig
	total_time = 50000.;
	N = (int)round(total_time/dt)+1;
	dt = (double)total_time/(N-1);  // correction dt
	
	dx = 0.025; //   0.01;
	L = 6.0;
	M = (int)round(L/dx)+1;
	dx = (double)L/(M-1);  // correction dx
	
	v = zero_vector(M);
	F = zero_vector(M);
	
	// points of measurement
	nb_times = 8;
	nb_space = 3; //15;
	st_times = zero_int_list(nb_times);
	st_indices = zero_int_list(nb_space);
	
	set_time_points(N);
	set_space_points(M);
	
	///////// parameters
	
	//fixed:
	//C_m = 1.;
	//kappa = 10.;
	pacing_width = 10;
	CV_bound = 10;
	
	//RD: ( Myocardial parameters)
	u_o = 0.0;
	u_u = 1.61;
	theta_v = 0.3;
	theta_w = 0.13;
	theta_v_m = 0.1;
	theta_o = 0.005;
	tau_v_1m = 80.0;
	tau_v_2m = 1.4506;
	tau_v_p = 1.4506;
	tau_w_1m = 70.0;
	tau_w_2m = 8.0;
	k_w = 200.0;
	u_w = 0.016;
	tau_w_p = 280.0;
	tau_fi = 0.078;
	tau_o1 = 410.0;
	tau_o2 = 7.0;
	tau_so1 = 91.0;
	tau_so2 = 0.8;
	k_so = 2.1;
	u_so = 0.6;
	tau_s1 = 2.7342;
	tau_s2 = 4.0;
	k_s = 2.0994;
	u_s = 0.9087;
	tau_si = 3.3849;
	tau_w_inf = 0.01;
	w_inf = 0.5;
	D = 0.001171;
	
	// calcium
	
	a_c = 1./7.;
	c0 = 0.01;
	tau_c = 60.;
	alpha = 0.004;
	c50 = 1.;
	rho0 = 0.1;
	rho1 = 0.163;
	Gamma = 2.6;
	bmax = 2.26;
	
	// mechanics
	
	U_s = 1.0;  // 0.95
	G_s = 0.1;//0.0013;
	F_T = 5.0;
	E = 1.;
	C1 = E/12.;
	C2 = E/12.;   // sigma = 6(c1+c2)*du/dx
	tau_T = 10.;
	tau_T2 = 1.;
	u_T = 0.05;
	
	E_o = E;
	l_o = L/5.;
	a = 1./(1.+2.*E*l_o/(E_o*L));
	
	// pacing
	pacing_Q = 1.;
	T_pace1 = 1000.0;
	T_pace2 = 1000.0;
	pacing_protocol = DYNAMIC;
	
	// measures
	thr_u = 0.2;
	eps_APD = 0.0001;
	eps_Newton = 1E-6;
	equilibrium_flag = 0;
	block_flag = 0;
	glob_flag = 0;
	S1S2_counter = 0;
	
	// others 
	finish_mode = FIXED;
}

void init_sim(){
	int i;
	
	T = (double**)malloc(M*sizeof(double*));
	DI = (double**)malloc(M*sizeof(double*));
	APD = (double**)malloc(M*sizeof(double*));
	CV = (double**)malloc(M*sizeof(double*));
	S2_index_DI = zero_int_list(M);
	S2_index_APD = zero_int_list(M);
	S2_index_CV = zero_int_list(M);
	S2_index_T = zero_int_list(M);
	Last_crossing = zero_vector(M);
	u = zero_vector(M);
	v = zero_vector(M);
	w = zero_vector(M);
	s = zero_vector(M);
	c = zero_vector(M);
	b = zero_vector(M);
	z = zero_vector(M);
	Ta = zero_vector(M);
	F = zero_vector(M);
	Old_u = zero_vector(M);
	for (i=0;i<M;i++){
		x = (double)i/(M-1)*L;
		u[i] = fu_init(x);
		v[i] = fv_init(x);
		w[i] = fw_init(x);
		s[i] = fs_init(x);
		c[i] = fc_init(x);
		b[i] = fb_init(x);
		z[i] = fz_init(x);
		Ta[i] = fT_init(x);
		F[i] = fF_init(x);
		T[i] = NULL;
		DI[i] = NULL;
		APD[i] = NULL;
		CV[i] = NULL;
		Last_crossing[i] = NAN;
	}
	t = 0.;
	if (pacing_protocol==S1S2_FINISH){
		pacing_protocol = S1S2;
		S1S2_counter = 0;
	}
	equilibrium_flag = 0;
	block_flag = 0;
	glob_flag = 0;
	S1S2_index_DI = 0;
	S1S2_index_APD = 0;
	T_size = zero_int_list(M);
	DI_size = zero_int_list(M);
	APD_size = zero_int_list(M);
	CV_size = zero_int_list(M);
	t_pace = T_pace1;
}

void finish(){
	int i;
	if (u!=NULL) free(u);
	if (v!=NULL) free(v);
	if (w!=NULL) free(w);
	if (s!=NULL) free(s);
	if (c!=NULL) free(c);
	if (b!=NULL) free(b);
	if (z!=NULL) free(z);
	if (F!=NULL) free(F);
	if (Ta!=NULL) free(Ta);
	if (Old_u!=NULL) free(Old_u);
	if (Last_crossing==NULL) free(Last_crossing);
	for (i=0;i<M;i++){
		if (T[i]!=NULL) free(T[i]);
		if (DI[i]!=NULL) free(DI[i]);
		if (APD[i]!=NULL) free(APD[i]);
		if (CV[i]!=NULL) free(CV[i]);
	}
	free(T);
	free(DI);
	free(APD);
	free(CV);
	free(T_size);
	free(DI_size);
	free(CV_size);
	free(APD_size);
	free(S2_index_DI);
	free(S2_index_APD);
	free(S2_index_CV);
	free(S2_index_T);
}

void init_data_storage(){
	char S[512];
	char ST[512];
	char SS[512];
	char Name[512];
	sprintf(S,"%s/%s",output_dir,output_name);
	if (mkdir(S,0777)==-1){
		printf("file %s exists allready -> abort...\n",S);
		exit(0);
	}
	if (nb_space>0){
		int i;
		sprintf(ST,"%s/%s/time",output_dir,output_name);
		mkdir(ST,0777);
		Time_series = (FILE**)malloc(nb_space*sizeof(FILE*));
		for (i=0;i<nb_space;i++){
			sprintf(Name,"%s/x%f",ST,(double)st_indices[i]/(M-1)*L);
			Time_series[i] = fopen(Name,"w");
		}
	}
	if (nb_times>0){
		sprintf(SS,"%s/%s/space",output_dir,output_name);
		mkdir(SS,0777);
	}
	sprintf(S,"%s/%s/restitution",output_dir,output_name);
	mkdir(S,0777);
	sprintf(S,"%s/%s/spacetime_u",output_dir,output_name);
	Spacetime_u = fopen(S,"w");
}

void finish_data_storage(){
	int i;
	for (i=0;i<nb_space;i++) if (Time_series[i]!=NULL) fclose(Time_series[i]);
	free(Time_series);
	Time_series = NULL;
	fclose(Spacetime_u);
}

void save_time(int k,int interval){
	int i,n;
	if ((k % interval)==0){
		for (i=0;i<nb_space;i++){
			n = st_indices[i];
			fprintf(Time_series[i],"%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\n",t,u[n],v[n],w[n],s[n],c[n],b[n],z[n],Ta[n],F[n]);
		}
	}
}

void save_space(int n,int interval,char* Name){
	int i;
	for (i=0;i<nb_times;i++) if (st_times[i]==n){
		FILE* file;
		char Full_name[512];
		int j;
		double x;
		if (Name==NULL) sprintf(Full_name,"%s/%s/space/t%f",output_dir,output_name,(double)n/(N-1)*total_time);
		else sprintf(Full_name,"%s/%s/%s",output_dir,output_name,Name);
		file = fopen(Full_name,"w");
		for (j=0;j<M;j++) if ((j % interval)==0){
			x = (double)j/(M-1)*L;
			fprintf(file,"%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\n",x,u[j],v[j],w[j],s[j],c[j],b[j],z[j],Ta[j],F[j]);
		}
		fclose(file);
	}
}

/*void save_initial_data(){
	int i;
	FILE* file;
	char Full_name[512];
	sprintf(Full_name,"%s/%s/initial",output_dir,output_name);
	file = fopen(Full_name,"w");
	fprintf(file,"%d\n",M);
	for (i=0;i<M;i++){
		fprintf(file,"%f\t%f\t%f\t%f\t%f\n",u[j],v[j],w[j],Ta[j],F[j]);
	}
	fclose(file);
}

void load_intial_data(char* Name){
	
}*/

void load_space_data(char* Dir,char* Name,char* Info_name){
	int ix;
	double X;
	char S[512];
	char* Col;
	char Line[1024];
	char Buffer[1024];
	char* read_success;
	sprintf(S,"%s/%s",Dir,Name);
	FILE* file = fopen(S,"r");
	if (file!=NULL){
		printf("load file %s ...\n",S);
		do{
			read_success = fgets(Line,1024,file);
			if (read_success!=NULL){
				sprintf(Buffer,"%s",Line);
				Col = strtok(Line,"\t\n");
				X = atof(Col);				
				ix = (int)round(X/dx);
				Col = strtok(NULL,"\t\n");
				u[ix] = atof(Col);	
				Col = strtok(NULL,"\t\n");
				v[ix] = atof(Col);	
				Col = strtok(NULL,"\t\n");
				w[ix] = atof(Col);		
				Col = strtok(NULL,"\t\n");
				s[ix] = atof(Col);		
				Col = strtok(NULL,"\t\n");
				c[ix] = atof(Col);		
				Col = strtok(NULL,"\t\n");
				b[ix] = atof(Col);		
				Col = strtok(NULL,"\t\n");
				z[ix] = atof(Col);
				Col = strtok(NULL,"\t\n");
				Ta[ix] = atof(Col);		
				Col = strtok(NULL,"\t\n");
				F[ix] = atof(Col);
			}
		}while(read_success!=NULL);
		fclose(file);
		sprintf(S,"%s/%s",output_dir,Info_name);
		file = fopen(S,"r");
		if (file!=NULL){		
			read_success = fgets(Line,1024,file);
			if (read_success!=NULL){
				sprintf(Buffer,"%s",Line);
				Col = strtok(Line,"=");				
				Col = strtok(NULL,"=");
				t = 0;//atof(Col);						Achtung: t=0 ms gesetzt !
			}	
			read_success = fgets(Line,1024,file);
			if (read_success!=NULL){
				sprintf(Buffer,"%s",Line);
				Col = strtok(Line,"=");				
				Col = strtok(NULL,"=");
				S1S2_counter = atoi(Col);
			}
			read_success = fgets(Line,1024,file);
			if (read_success!=NULL){
				sprintf(Buffer,"%s",Line);
				Col = strtok(Line,"=");				
				Col = strtok(NULL,"=");
				t_pace = atof(Col);
			}						
			printf("completed\n");
			fclose(file);
		}
	}
	else{
		printf("error opening file %s -> abort\n",S);
		exit(0);
	}
}

void save_spacetime(int k){
	if ((k % (2*data_interval))==0){
		int j;
		char* Line;
		char Col[512];
		Line = (char*)malloc(32*M*sizeof(char));
		Line[0]='\0';
		for (j=0;j<M-1;j++){
			sprintf(Col,"%f\t",u[j]);
			strcat(Line,Col);
		}
		sprintf(Col,"%f\n",u[M-1]);
		strcat(Line,Col);
		fprintf(Spacetime_u,"%s",Line);
		free(Line);
	}
}

int APD_equilibrium(int min){
	int i;
	double q;
	int flag = 1;
	for (i=0;i<M;i++){
		if (APD_size[i]>min){
			q = (APD[i][APD_size[i]-1]-APD[i][APD_size[i]-2])/APD[i][APD_size[i]-1];
			if (fabs(q)>eps_APD) flag = 0;
		}
		else flag = 0;
	}
	return flag;
}

int APD_2nd_equilibrium(int min){
	int i;
	double q1,q2,max1,max2;
	int flag = 1;
	max1 = 0;
	max2 = 0;
	for (i=0;i<M;i++){
		if (APD_size[i]>min){
			q1 = (APD[i][APD_size[i]-1]-APD[i][APD_size[i]-3])/APD[i][APD_size[i]-1];
			q2 = (APD[i][APD_size[i]-2]-APD[i][APD_size[i]-4])/APD[i][APD_size[i]-2];
			if (q1>max1) max1 = q1;
			if (q2>max2) max2 = q2;
			if (fabs(q1)>eps_APD || fabs(q2)>eps_APD) flag = 0;
		}
		else flag = 0;
	}
	//printf("\neps1=%e max2=%e\n",max1,max2); 
	return flag;
}

void S1S2_identification(){
	int i,j;
	double d,s,di,apd,cv,tp,fac,X;
	char Name[512];
	sprintf(Name,"%s/%s/restitution/S1S2",output_dir,output_name);
	FILE* file = fopen(Name,"w");
	if (file!=NULL){
		for (i=CV_bound;i<M;i++){
			fac = 50.;//1./eps_APD;
			X = (double)i/(M-1)*L;
			fprintf(file,"%f\t",X);
			di = NAN;
			d = fabs(DI[i][S2_index_DI[i]]-DI[i][S2_index_DI[i]-1]);
			for (j=S2_index_DI[i]+1;j<DI_size[i];j++){
				s = fabs(DI[i][j]-DI[i][j-1]);
				if ((s/d)>fac){
					di = DI[i][j];
					break;
				}
			}
			fprintf(file,"%f\t",di);
			apd = NAN;
			cv = NAN;
			tp = NAN;
			d = fabs(APD[i][S2_index_APD[i]]-APD[i][S2_index_APD[i]-1]);
			for (j=S2_index_APD[i]+1;j<APD_size[i];j++){
				s = fabs(APD[i][j]-APD[i][j-1]);
				if ((s/d)>fac){
					apd = APD[i][j];
					cv = CV[i][j];
					tp = T[i][j];
					break;
				}
			}
			fprintf(file,"%f\t%f\t%f\n",apd,cv,tp);				
		}
		fclose(file);
	}
	else{
		printf("error creating file %s -> abort",Name);
		exit(0);
	}
}

void S1S2_save(){
	int i;
	double di,apd,cv,tp,X;
	char Name[512];
	sprintf(Name,"%s/%s/restitution/S1S2",output_dir,output_name);
	FILE* file = fopen(Name,"w");
	if (file!=NULL){
		for (i=CV_bound;i<M;i++){			
			X = (double)i/(M-1)*L;			
			apd = APD[i][S1S2_index_APD+1];
			di = DI[i][S1S2_index_DI+1];
			cv = CV[i][S1S2_index_DI+1];
			tp = T[i][S1S2_index_DI+1];			
			fprintf(file,"%f\t%f\t%f\t%f\t%f\n",X,di,apd,cv,tp);				
		}
		fclose(file);
	}
	else{
		printf("error creating file %s -> abort",Name);
		exit(0);
	}
}

void dynamic_save(){
	int i;
	double X;
	char Name[512];
	sprintf(Name,"%s/%s/restitution/DYN",output_dir,output_name);
	FILE* file = fopen(Name,"w");
	int alternans = 1-APD_equilibrium(4);
	for (i=CV_bound;i<M;i++){
		X = (double)i/(M-1)*L;
		if (DI_size[i]>1 && APD_size[i]>1 && CV_size[i]>1 && T_size[i]>1){
			fprintf(file,"%f\t%f\t%f\t%f\t%f\n",X,DI[i][APD_size[i]-2],
			 APD[i][APD_size[i]-1],CV[i][CV_size[i]-1],T[i][T_size[i]-1]);
		}
		if (alternans && DI_size[i]>2 && APD_size[i]>2 && CV_size[i]>2 && T_size[i]>2){
			fprintf(file,"%f\t%f\t%f\t%f\t%f\n",X,DI[i][APD_size[i]-3],
		 	 APD[i][APD_size[i]-2],CV[i][CV_size[i]-2],T[i][T_size[i]-2]);
		}
	}
	fclose(file);
}

void save_rest_prop(){
	int i,j;
	double X;
	FILE* file;
	char Name[512];
	char Col[512];
	char Line[65536];
	sprintf(Name,"%s/%s/restitution/DIx",output_dir,output_name);
	file = fopen(Name,"w");
	for (i=0;i<M;i++) if (DI_size[i]>0){
		X = (double)i/(M-1)*L;
		sprintf(Line,"%f\t",X);
		for (j=0;j<DI_size[i]-1;j++){
			sprintf(Col,"%f\t",DI[i][j]);
			strcat(Line,Col);
		}
		sprintf(Col,"%f\n",DI[i][DI_size[i]-1]);
		strcat(Line,Col);
		fprintf(file,"%s",Line);
	}
	fclose(file);
	
	sprintf(Name,"%s/%s/restitution/APDx",output_dir,output_name);
	file = fopen(Name,"w");
	for (i=0;i<M;i++) if (APD_size[i]>0){
		X = (double)i/(M-1)*L;
		sprintf(Line,"%f\t",X);
		for (j=0;j<APD_size[i]-1;j++){
			sprintf(Col,"%f\t",APD[i][j]);
			strcat(Line,Col);
		}
		sprintf(Col,"%f\n",APD[i][APD_size[i]-1]);
		strcat(Line,Col);
		fprintf(file,"%s",Line);
	}
	fclose(file);
	
	sprintf(Name,"%s/%s/restitution/CVx",output_dir,output_name);
	file = fopen(Name,"w");
	for (i=CV_bound;i<M-CV_bound;i++) if (CV_size[i]>0){
		X = (double)i/(M-1)*L;
		sprintf(Line,"%f\t",X);
		for (j=0;j<CV_size[i]-1;j++){
			sprintf(Col,"%f\t",CV[i][j]);
			strcat(Line,Col);
		}
		sprintf(Col,"%f\n",CV[i][CV_size[i]-1]);
		strcat(Line,Col);
		fprintf(file,"%s",Line);
	}
	fclose(file);
	
	sprintf(Name,"%s/%s/restitution/Tx",output_dir,output_name);
	file = fopen(Name,"w");
	for (i=0;i<M;i++) if (T_size[i]>0){
		X = (double)i/(M-1)*L;
		sprintf(Line,"%f\t",X);
		for (j=0;j<T_size[i]-1;j++){
			sprintf(Col,"%f\t",T[i][j]);
			strcat(Line,Col);
		}
		sprintf(Col,"%f\n",T[i][T_size[i]-1]);
		strcat(Line,Col);
		fprintf(file,"%s",Line);
	}
	fclose(file);
	
	switch(pacing_protocol){
		case S1S2_FINISH: 			
			S1S2_save();
			break;
		case DYNAMIC:		
			dynamic_save();
			break;
	}
}

int take_all(const struct dirent* Entry){
	if (strcmp(Entry->d_name,"..")==0 || strcmp(Entry->d_name,".")==0) return 0;
	else return 1;
}

void save_rest_curves(char* Dir){
	int i,k,ix;
	double X;
	char* read_success;
	char* Col;
	char* Mode;
	char Name[512];
	char Line[512];
	char Buffer[512];
	FILE* Data;
	struct dirent** Entries;
	if (pacing_protocol==DYNAMIC) Mode = "DYN"; else Mode = "S1S2";
	int number = scandir(Dir,&Entries,take_all,alphasort);
	if (number>=0){
		FILE** Rest = (FILE**)malloc(nb_space*sizeof(FILE*));
		for (i=0;i<nb_space;i++){
			X = (double)st_indices[i]/(M-1)*L;
			sprintf(Name,"%s/rest_curve_x%f",Dir,X);
			Rest[i] = fopen(Name,"w");
		}
		for (i=0;i<number;i++){
			sprintf(Name,"%s/%s/restitution/%s",Dir,Entries[i]->d_name,Mode);
			Data = fopen(Name,"r");
			if (Data!=NULL){
				printf("open %s\n",Name);
				do{
					read_success = fgets(Line,512,Data);
					sprintf(Buffer,"%s",Line);
					if (read_success!=NULL){
						Col = strtok(Line,"\t");
						X = atof(Col);
						ix = (int)round(X/dx);
						for (k=0;k<nb_space;k++) if (ix==st_indices[k]){
							fprintf(Rest[k],"%s",Buffer);
						}
					}
				}while(read_success!=NULL);
				fclose(Data);
			}			
		}
		for (i=0;i<nb_space;i++) fclose(Rest[i]);
	}
}

void save_param(){
	char Name[512];
	sprintf(Name,"%s/%s/parameters",output_dir,output_name);
	FILE* file = fopen(Name,"w");
	fprintf(file,"dx=%f\n",dx);
	fprintf(file,"L=%f\n",L);
	fprintf(file,"M=%d\n",M);
	fprintf(file,"dt=%f\n",dt);
	if (pacing_protocol==DYNAMIC){
		fprintf(file,"protocol=dynamic\n");
		fprintf(file,"Tpace=%f\n",T_pace1);
	}
	else{
		fprintf(file,"protocol=S1S2\n");
		fprintf(file,"S1=%f\n",T_pace1);
		fprintf(file,"S2=%f\n",T_pace2);
	}
	fprintf(file,"u_o=%f\n",u_o);
	fprintf(file,"u_u=%f\n",u_u);
	fprintf(file,"theta_v=%f\n",theta_v);
	fprintf(file,"theta_w=%f\n",theta_w);
	fprintf(file,"theta_v_m=%f\n",theta_v_m);
	fprintf(file,"theta_o=%f\n",theta_o);
	fprintf(file,"tau_v_1m=%f\n",tau_v_1m);
	fprintf(file,"tau_v_2m=%f\n",tau_v_2m);
	fprintf(file,"tau_v_p=%f\n",tau_v_p);
	fprintf(file,"tau_w_1m=%f\n",tau_w_1m);
	fprintf(file,"tau_w_2m=%f\n",tau_w_2m);
	fprintf(file,"tau_w_p=%f\n",tau_w_p);
	fprintf(file,"k_w=%f\n",k_w);
	fprintf(file,"u_w=%f\n",u_w);
	fprintf(file,"tau_fi=%f\n",tau_fi);
	fprintf(file,"tau_o1=%f\n",tau_o1);
	fprintf(file,"tau_o2=%f\n",tau_o2);
	fprintf(file,"tau_so1=%f\n",tau_so1);
	fprintf(file,"tau_so2=%f\n",tau_so2);
	fprintf(file,"k_so=%f\n",k_so);
	fprintf(file,"u_so=%f\n",u_so);
	fprintf(file,"tau_s1=%f\n",tau_s1);
	fprintf(file,"tau_s2=%f\n",tau_s2);
	fprintf(file,"k_s=%f\n",k_s);
	fprintf(file,"u_s=%f\n",u_s);
	fprintf(file,"w_inf=%f\n",w_inf);
	fprintf(file,"tau_si=%f\n",tau_si);
	fprintf(file,"tau_w_inf=%f\n",tau_w_inf);
	
	fprintf(file,"a_c=%f\n",a_c);
	fprintf(file,"c0=%f\n",c0);
	fprintf(file,"tau_c=%f\n",tau_c);
	fprintf(file,"alpha=%f\n",alpha);
	fprintf(file,"c50=%f\n",c50);
	fprintf(file,"rho0=%f\n",rho0);
	fprintf(file,"rho1=%f\n",rho1);
	fprintf(file,"Gamma=%f\n",Gamma);
	fprintf(file,"bmax=%f\n",bmax);
	
	fprintf(file,"D=%f\n",D);
	fprintf(file,"U_s=%f\n",U_s);
	fprintf(file,"G_s=%f\n",G_s);
	fprintf(file,"F_T=%f\n",F_T);
	fprintf(file,"E=%f\n",E);
	fprintf(file,"tau_T=%f\n",tau_T);
	fprintf(file,"tau_T2=%f\n",tau_T2);
	fprintf(file,"u_T=%f\n",u_T);
	fprintf(file,"E_o=%f\n",E_o);
	fprintf(file,"l_o=%f\n",l_o);
	
	fclose(file);
}

void save_info(){
	char Name[512];
	sprintf(Name,"%s/%s/info",output_dir,output_name);
	FILE* file = fopen(Name,"w");
	fprintf(file,"time=%f\n",t);
	fprintf(file,"S1S2counter=%d\n",S1S2_counter);
	fprintf(file,"pacecountdown=%f\n",t_pace);
	fprintf(file,"equilibrium=%d\n",equilibrium_flag);
	fprintf(file,"block=%d\n",block_flag);
	fprintf(file,"globexcitation=%d\n",glob_flag);
	fclose(file);
}

double get_2nd_order_time(double y1,double y2,double y3,double y,double dt){  // Achtung: Nur wenn der Old_u auch geladen wurde sinnvoll !
	double cross;
	double rad = y1*y1+(-4.*y2+y3)*(-4.*y2+y3)+8.*y*(y1-2.*y2+y3)-2.*y1*(4.*y2+y3);
	if (rad>=0){
		cross = (y1-y3+sqrt(rad))*dt/(2.*(y1-2.*y2+y3));
		if (cross>0 && cross<dt) return cross;
		cross = (y1-y3-sqrt(rad))*dt/(2.*(y1-2.*y2+y3));
		if (cross>0 && cross<dt) return cross;
		return NAN;
	}
	else{
		printf("error in interpolation -> abort\n");
		exit(0);
	}
}

void set_rest_prop(double* U_new,double* U,double* U_old,double t){
	int i;
	for (i=0;i<M;i++) if ((U_new[i]-thr_u)*(U[i]-thr_u)<0){
		double t_cross = t+(thr_u-U[i])/(U_new[i]-U[i])*dt;					 //linear interpolation
		//double t_cross = t+get_2nd_order_time(U_old[i],U[i],U_new[i],thr_u,dt);  // quadratic interpolation
		if (!isnan(Last_crossing[i])){
			if (U_new[i]>U[i]){
				DI_size[i]++;
				DI[i] = (double*)realloc(DI[i],DI_size[i]*sizeof(double));
				DI[i][DI_size[i]-1] = t_cross-Last_crossing[i];
				if (APD_size[i]>0 && DI_size[i]>0){
					T_size[i]++;
					T[i] = (double*)realloc(T[i],T_size[i]*sizeof(double));
					T[i][T_size[i]-1] = APD[i][APD_size[i]-1]+DI[i][DI_size[i]-1];
				}
				if (i>CV_bound-1){
					double cv;
					cv = dx/(t_cross-Last_crossing[i-1]);
					CV_size[i]++;
					CV[i] = (double*)realloc(CV[i],CV_size[i]*sizeof(double));
					CV[i][CV_size[i]-1] = cv;
				}
			}
			else{				
				APD_size[i]++;
				APD[i] = (double*)realloc(APD[i],APD_size[i]*sizeof(double));
				APD[i][APD_size[i]-1] = t_cross-Last_crossing[i];								
			}
		}
		Last_crossing[i] = t_cross;
	}
}

int glob_excitation(double cmax,double t,double t0){		//-1: no excitation for some i
	int i;													// 0: normal conduction
	double Dt,cm;											// 1: global excitation
	double sum = 0;
	double X = 0;
	i = 0;
	if (t>t0){
		for (i=0;i<M;i++) if (isnan(Last_crossing[i])) return -1;
		for(i=1;i<M;i++){
			Dt = Last_crossing[i]-Last_crossing[i-1];
			sum += Dt;
			X += dx;
			if (Dt<0){
				cm = X/sum;
				if (u[i]>thr_u && (sum==0 || cm>cmax)){
					//printf("cmean=%f at t=%f\n",cm,t);
					//for (j=0;j<M;j++) printf("%f\n",Last_crossing[j]);
					return 1;					
				}
				sum = 0;
				X = 0;
			}
		}
		cm = X/sum;
		if (u[i]>thr_u && (sum==0 || cm>cmax)){
			//printf("cmean=%f at t=%f\n",cm,t);
			//for (j=0;j<M;j++) printf("%f\n",Last_crossing[j]);
			return 1;		
		}
		else return 0;
	}
	else return 0;
}

int block_detected(double t,double tmax){
	int i;
	for (i=0;i<M;i++){
		if (!isnan(Last_crossing[i])){
			if ((t-Last_crossing[i])>tmax) return 1;
		}
		else{
			if (t>tmax) return 1;	
		}		
	}
	return 0;
}

double mean_val(double* X,int n){
	int i;
	double res = 0;
	for (i=0;i<n;i++) res += X[i];
	return (double)res/n;
}

double my_pow(double x,int n){
	int i;
	double res = 1;
	for (i=0;i<n;i++) res *= x;
	return res;
}

double arb_power(double x, double a){
	return exp(a*log(x));
}

double tf(double x){
	return sin(x);
}

double sekant_root(double (*F)(double x),double x0,double x1,double eps){
	double z_new;
	int max = 100;
	double z_old = x0;
	double z = x1;
	double f0 = (*F)(z);
	int i = 0;
	do{
		z_new = z-(*F)(z)*(z-z_old)/((*F)(z)-(*F)(z_old));
		z_old = z;
		z = z_new;
		i++; 
	}while(fabs((*F)(z)/f0)>eps && i<max);
	if (i==max) printf("warning max iteration number reached ...\n");
	return z;
}

double integral(double (*F)(double x)){
	int i;
	double X;
	double sum = 0;
	for (i=0;i<M;i++){
		X = (double)i/(M-1)*L;
		sum += (*F)(X)*dx;
	}
	return sum;
}

double total_P(double F){
	return 2.*(C1*(F*F-1./F)+C2*(F-1./(F*F))+P_Ta/2.)-P_a*F;
}

double getF(double Ta,double a,double F0){
	double dF = 0.1;
	P_Ta = Ta;
	P_a = a;
	return sekant_root(&total_P,F0-dF,F0,eps_Newton);
}

double integral_F(double a){
	int i;
	//double X;
	double sum = 0;
	for (i=0;i<M;i++){
		//X = (double)i/(M-1)*L;
		sum += getF(Ta[i],a,F[i])*dx;
	}
	return (sum/L-1.);
}

double geta(double last_a){
	double da = 0.1;
	return sekant_root(&integral_F,last_a-da,last_a,eps_Newton);
}

double Heavyside(double x,double a){
	if (x<a) return 0; else return 1.;
}

double tau_o(double u){
	if (u<theta_o) return tau_o1; else return tau_o2;
}

double tau_s(double u){
	if (u<theta_w) return tau_s1; else return tau_s2;
}

double tau_so(double u){
	return tau_so1+(tau_so2-tau_so1)*(1.+tanh(k_so*(u-u_so)))/2.;
}

double tau_v_m(double u){
	if (u<theta_v_m) return tau_v_1m; else return tau_v_2m;
}

double v_inf(double u){
	if (u<theta_v_m) return 1.; else return 0.;
}

double fw_inf(double u){
	if (u<theta_o) return 1.-u/tau_w_inf; else return w_inf;
}

double tau_w_m(double u){
	return tau_w_1m+(tau_w_2m-tau_w_1m)*(1.+tanh(k_w*(u-u_w)))/2.;
}

double J_fi(double u,double v){
	return Heavyside(u,theta_v)*v*(u-u_u)*(u-theta_v)/tau_fi;
}

double J_si(double u,double w,double s){
	return -Heavyside(u,theta_w)*w*s/tau_si;
}

double J_so(double u){
	return (u-u_o)*(1.-Heavyside(u,theta_w))/tau_o(u)+Heavyside(u,theta_w)/tau_so(u);
}

double J_sa(double u,double T,double T_mean){
	return G_s*(u-U_s)*(a*T_mean-T)*Heavyside(a*T_mean,T)/E;
}

//double f_tau(double u){
	//return tau_v2+(tau_v1-tau_v2)*Heavyside(u,u_v);
//}

double f_U(double u,double v,double w, double s){
	return -J_fi(u,v)-J_so(u)-J_si(u,w,s);
}

double f_V(double u,double v){
	return -v/tau_v_p+Heavyside(theta_v,u)*((v_inf(u)-v)/tau_v_m(u)+v/tau_v_p);
}

double f_W(double u,double w){
	return -w/tau_w_p+Heavyside(theta_w,u)*((fw_inf(u)-w)/tau_w_m(u)+w/tau_w_p);
}

double f_S(double u,double s){
	return ((1.+tanh(k_s*(u-u_s)))/2.-s)/tau_s(u);
}

double f_C(double u,double w, double s, double c){
	return -a_c*J_si(u,w,s)-(c-c0)/tau_c;
}

double f_B(double c, double b){
	return rho0*c*(bmax-b)-rho1*(1.-1./Gamma)*b;
}

double f_Z(double b, double z){
	return alpha*(arb_power(b/c50,4.5)*(1.-z)-z);
}

double get_T(double z){
	return F_T*z;
}

double f_T(double u,double T){
	if (u>u_T) return (F_T*u-T)/tau_T;
	else return (F_T*u-T)/tau_T2;
}

/*double f_T(double u, double T){
	double APD_max = 450.;
	double b = 50.;
	double a = F_T/APD_max;
	if (u>u_T) return a*u;
	else return (F_T*u-T)/b;
}*/

/*double f_T(double u,double w, double s, double T){
	tau_T = 50.0;
	return -J_si(u,w,s)-T/tau_T;
}*/

double f_sa(double F,double u){   // dummy
	return 0;
}

void pacing(double* U,int finish_mode){
	if (t_pace>=T_pace1){
		int i;
		for (i=0;i<pacing_width;i++) U[i] += pacing_Q;
		t_pace = 0.;
		
		/*if (block_detected(t,10.*T_pace1)){
			block_flag = 1;
			printf("conduction block detected -> finish\n");
		}
		
		if (glob_excitation(10.,t,10.*T_pace1)==1){
			glob_flag = 1;
			printf("global ecxitation detected -> finish\n");			
		}*/
		
		if (pacing_protocol==S1S2_FINISH){
			if (S1S2_counter==3){
				equilibrium_flag =1;
				printf("\n S1S2 finished\n");
			}
			S1S2_counter++;
		}
		
		if (APD_equilibrium(4) || APD_2nd_equilibrium(4)){
			if (pacing_protocol==S1S2){
				printf("\n equilibrium reached\n");
				printf("set S2 pacing\n");
				T_pace1 = T_pace2;
				S1S2_index_DI = DI_size[0];
				S1S2_index_APD = APD_size[0];
				for (i=1;i<M;i++){
					S2_index_DI[i] = DI_size[i]-1;
					S2_index_APD[i] = APD_size[i]-1;
					S2_index_CV[i] = CV_size[i]-1;
					S2_index_T[i] = T_size[i]-1;
				}
				pacing_protocol = S1S2_FINISH;
			}
			if (pacing_protocol==DYNAMIC && finish_mode==EQUILIBRIUM){
				equilibrium_flag = 1;
				printf("\n equilibrium reached -> finish\n");
			}
		}
	}
	else t_pace += dt;
}



void next_time_step(double* U,double* V,double* W,double* S,double* C, double* B, double* Z, double* TA,double* F,double t){ // explicit Euler 
	int i;
	double q = D*dt/(dx*dx);
	double* New_u = clone_vector(U,M);
	double mean_T = mean_val(TA,M);
	//static double a = 0;
	
	// diffusion for u
	New_u[0] += q*(U[1]-U[0]);  // Neumann BC for u
	for (i=1;i<M-1;i++) New_u[i] += q*(U[i+1]-2.*U[i]+U[i-1]);
	New_u[M-1] += q*(U[M-2]-U[M-1]);
	
	// reaction 
	for (i=0;i<M;i++){
		New_u[i] += (f_U(U[i],V[i],W[i],S[i])-J_sa(U[i],TA[i],mean_T))*dt;
		V[i] += f_V(U[i],V[i])*dt;
		W[i] += f_W(U[i],W[i])*dt;
		S[i] += f_S(U[i],S[i])*dt;
		C[i] += f_C(U[i],W[i],S[i],C[i])*dt;
		B[i] += f_B(C[i],B[i])*dt;
		Z[i] += f_Z(B[i],Z[i])*dt;
		TA[i] = get_T(Z[i]);//f_T(U[i],TA[i])*dt;;
		F[i] = (a*mean_T-TA[i])/E;//J_sa(U[i],TA[i],mean_T);
	}
	
	// mechanics
	/*a = geta(a);
	printf("a=%f\n",a);*/
	/*for (i=0;i<M;i++){
		F[i] = 1.;//getF(TA[i],a,F[i]);
	}*/
	
	pacing(New_u,finish_mode);
	set_rest_prop(New_u,U,Old_u,t);
	for (i=0;i<M;i++){
		Old_u[i] = U[i];
		U[i] = New_u[i];
	}
	free(New_u);
}

int simulate(char* load,char* store){
	int m;
	int i = 0;
	t = 0;
	init_sim();
	if (load!=NULL){
		char Name[512];
		char Info[512];
		sprintf(Name,"%s/initial",load);
		sprintf(Info,"%s/info",load);
		load_space_data(output_dir,Name,Info);
	}
	init_data_storage();
	save_param();
	
	// simulation
	printf("start simulation: output -> %s\n",output_name);
	
	while(i<N && !equilibrium_flag && !block_flag && !glob_flag){
		save_time(i,data_interval);
		save_space(i,space_interval,NULL);
		if (i>(9*N/10)) save_spacetime(i);
		next_time_step(u,v,w,s,c,b,z,Ta,F,t);
		if (i>0 && i % (N/100)== 0){
			m = 100*i/N;
			printf("\r%d%% abgeschlossen",m);
			fflush(stdout);
		}
		if (store!=NULL && pacing_protocol==S1S2_FINISH) equilibrium_flag = 1;  // force exit for S1S2 in store mode
		t += dt;
		i++;
	}
	i--;
	t -=dt;
	
	// end simulation
	if (store!=NULL) save_space(i,1,store);
	save_rest_prop();
	save_info();
	finish_data_storage();
	finish();
	printf("\n finished !\n");
	return 0;
}

void compute_rest_curve(double Tmin,double Tmax,double S1,int n){
	int i;
	char Dir[512];
	char Orig[512];
	//char Name[512];
	sprintf(Orig,"%s",output_name);
	sprintf(Dir,"%s/%s",output_dir,output_name);
	if (mkdir(Dir,0777)==-1){
		printf("file %s exists allready -> abort...\n",Dir);
		exit(0);
	}
	sprintf(output_dir,"%s",Dir);
	finish_mode = EQUILIBRIUM;
	for (i=0;i<n;i++){
		sprintf(output_name,"%s_%d",Orig,i);
		if (pacing_protocol==DYNAMIC) T_pace1 = (double)Tmin+(Tmax-Tmin)*i/(n-1);
		if (pacing_protocol==S1S2 || pacing_protocol==S1S2_FINISH){
			T_pace1 = S1;
			T_pace2 = (double)Tmin+(Tmax-Tmin)*i/(n-1);
		}
		simulate(NULL,NULL);
	}
}

/*void compute_rest_curve(double Tmin,double Tmax,double S1,int n){
	int i;
	char Dir[512];
	char Orig[512];
	//char Name[512];
	sprintf(Orig,"%s",output_name);
	sprintf(Dir,"%s/%s",output_dir,output_name);
	if (mkdir(Dir,0777)==-1){
		printf("file %s exists allready -> abort...\n",Dir);
		exit(0);
	}
	sprintf(output_dir,"%s",Dir);
	finish_mode = EQUILIBRIUM;
	if (pacing_protocol==DYNAMIC){
		for (i=0;i<n;i++){
			sprintf(output_name,"%s_%d",Orig,i);
			T_pace1 = (double)Tmin+(Tmax-Tmin)*i/(n-1);
			simulate(NULL,NULL);
		}
	}
	else{
		T_pace1 = S1;
		simulate(NULL,"last");
		for (i=0;i<n;i++){
			sprintf(output_name,"%s_%d",Orig,i);
			T_pace1 = S1;
			T_pace2 = (double)Tmin+(Tmax-Tmin)*i/(n-1);
			simulate("last",NULL);
		}
	}
}*/

void vary_gs(double gmin,double gmax,int n,double Length,double T_pace,int continuation){
	int i;
	char Dir[512];
	char Orig[512];
	char* Prev;
	char* Next;
	//char Name[512];
	sprintf(Orig,"%s",output_name);
	sprintf(Dir,"%s/%s",output_dir,output_name);
	if (mkdir(Dir,0777)==-1){
		printf("file %s exists allready -> abort...\n",Dir);
		exit(0);
	}
	sprintf(output_dir,"%s",Dir);
	finish_mode = FIXED;
	pacing_protocol = DYNAMIC;
	T_pace1 = T_pace;
	L = Length;
	M = (int)round(L/dx)+1;
	dx = (double)L/(M-1);  // correction dx	
	set_space_points(M);
	for (i=0;i<n;i++){
		sprintf(output_name,"%s_%d",Orig,i);
		G_s = (double)gmin+(gmax-gmin)*i/(n-1);
		if (continuation==INDEPENDENT || i==0) Prev = NULL; 
		else{
			Prev = (char*)malloc(512*sizeof(char));
			sprintf(Prev,"%s_%d",Orig,i-1);
		}
		if (continuation==INDEPENDENT || i==n-1) Next = NULL; 
		else{
			Next = "initial";
		}
		simulate(Prev,Next);
	}
}

int main(int argc, char* argv[]){
	
	
	/*int i;
	double X;
	int K = 65536;
	complex* Data = (complex*)malloc(K*sizeof(complex));
	for (i=0;i<K;i++){
		X = (double)i/(K-1)*1.;
		Data[i].x = X*(X-1.);
		Data[i].y = 0;
	}
	init_FastFourier1D(Data,K);
	set_FFT_mode(ANALYSE);
	fastfourier();
	complex* Trans = cmp_cpy(get_FFT_trans(),K);
	init_FastFourier1D(Trans,K);
	set_FFT_mode(SYNTHESE);
	fastfourier();
	complex* N_Trans = get_FFT_trans();
	
	exit(0);*/
	
	/*double x0 = 4.;
	double x1 = 3.5;
	double z = sekant_root(&tf,x0,x1,0.0001);
	printf("root: %f\n",z);
	exit(0);*/
	

	/*int i;	
	double X;
	M = 100;
	L = 1.;
	dx = (double)L/M; 
	Ta = zero_vector(M);
	F = zero_vector(M);
	for (i=0;i<M;i++){
		X = (double)i/(M-1)*L;
		Ta[i] = exp(-10.*(x-0.5)*(x-0.5))/5.;
		F[i] = 1.;
	}
	double a0 = 1.;
	double a = geta(a0);
	exit(0);*/
	
	char Full_name[512];
	// init data aquisition
	sprintf(output_dir,"/users/radszuweit/Daten");
	sprintf(output_name,"%s",argv[1]);
	sprintf(Full_name,"%s/%s",output_dir,output_name);
	double Length;
	double period;
	init_param();
	if (argc>2){
		Length = atof(argv[2]);
		period = atof(argv[3]);
		//vary_gs(0.000,0.005,20,Length,period,CONTINUATION);
		vary_gs(0.0,0.006,15,Length,period,CONTINUATION);
	}
	
	// main program
	//compute_rest_curve(350,1500,1500,30);
	simulate(NULL,NULL);
	//save_rest_curves(Full_name);
	
	// clean
	free(st_times);
	free(st_indices);
	
	return 0;
}

