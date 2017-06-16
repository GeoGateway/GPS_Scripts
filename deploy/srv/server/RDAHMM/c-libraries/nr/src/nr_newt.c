#include <math.h>
#define NRANSI
#include "nr_util.h"
#define MAXITS 200
#define TOLF 1.0e-4
#define TONR_LMIN 1.0e-6
#define TOLX 1.0e-7
#define STPMX 100.0

int nn;
float *fvec;
void (*nrfuncv)(int n, float v[], float f[]);
#define FREERETURN {NR_free_vector(fvec,1,n);NR_free_vector(xold,1,n);\
	NR_free_vector(p,1,n);NR_free_vector(g,1,n);NR_free_matrix(fjac,1,n,1,n);\
	NR_free_ivector(indx,1,n);return;}

void NR_newt(float x[], int n, int *check,
	void (*vecfunc)(int, float [], float []))
{
	void NR_fdjac(int n, float x[], float fvec[], float **df,
		void (*vecfunc)(int, float [], float []));
	float NR_fmin(float x[]);
	void NR_lnsrch(int n, float xold[], float fold, float g[], float p[], float x[],
		 float *f, float stpmax, int *check, float (*func)(float []));
	void NR_lubksb(float **a, int n, int *indx, float b[]);
	void NR_ludcmp(float **a, int n, int *indx, float *d);
	int i,its,j,*indx;
	float d,den,f,fold,stpmax,sum,temp,test,**fjac,*g,*p,*xold;

	indx=NR_ivector(1,n);
	fjac=NR_matrix(1,n,1,n);
	g=NR_vector(1,n);
	p=NR_vector(1,n);
	xold=NR_vector(1,n);
	fvec=NR_vector(1,n);
	nn=n;
	nrfuncv=vecfunc;
	f=NR_fmin(x);
	test=0.0;
	for (i=1;i<=n;i++)
		if (fabs(fvec[i]) > test) test=fabs(fvec[i]);
	if (test < 0.01*TOLF) {
		*check=0;
		FREERETURN
	}
	for (sum=0.0,i=1;i<=n;i++) sum += NR_SQR(x[i]);
	stpmax=STPMX*NR_FMAX(sqrt(sum),(float)n);
	for (its=1;its<=MAXITS;its++) {
		NR_fdjac(n,x,fvec,fjac,vecfunc);
		for (i=1;i<=n;i++) {
			for (sum=0.0,j=1;j<=n;j++) sum += fjac[j][i]*fvec[j];
			g[i]=sum;
		}
		for (i=1;i<=n;i++) xold[i]=x[i];
		fold=f;
		for (i=1;i<=n;i++) p[i] = -fvec[i];
		NR_ludcmp(fjac,n,indx,&d);
		NR_lubksb(fjac,n,indx,p);
		NR_lnsrch(n,xold,fold,g,p,x,&f,stpmax,check,NR_fmin);
		test=0.0;
		for (i=1;i<=n;i++)
			if (fabs(fvec[i]) > test) test=fabs(fvec[i]);
		if (test < TOLF) {
			*check=0;
			FREERETURN
		}
		if (*check) {
			test=0.0;
			den=NR_FMAX(f,0.5*n);
			for (i=1;i<=n;i++) {
				temp=fabs(g[i])*NR_FMAX(fabs(x[i]),1.0)/den;
				if (temp > test) test=temp;
			}
			*check=(test < TONR_LMIN ? 1 : 0);
			FREERETURN
		}
		test=0.0;
		for (i=1;i<=n;i++) {
			temp=(fabs(x[i]-xold[i]))/NR_FMAX(fabs(x[i]),1.0);
			if (temp > test) test=temp;
		}
		if (test < TOLX) FREERETURN
	}
	NR_error("MAXITS exceeded in NR_newt");
}
#undef MAXITS
#undef TOLF
#undef TONR_LMIN
#undef TOLX
#undef STPMX
#undef FREERETURN
#undef NRANSI
