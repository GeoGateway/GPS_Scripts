#include <math.h>
#define NRANSI
#include "nr_util.h"

void usrfun(float *x,int n,float *fvec,float **fjac);
#define FREERETURN {NR_free_matrix(fjac,1,n,1,n);NR_free_vector(fvec,1,n);\
	NR_free_vector(p,1,n);NR_free_ivector(indx,1,n);return;}

void NR_mNR_newt(int ntrial, float x[], int n, float tolx, float tolf)
{
	void NR_lubksb(float **a, int n, int *indx, float b[]);
	void NR_ludcmp(float **a, int n, int *indx, float *d);
	int k,i,*indx;
	float errx,errf,d,*fvec,**fjac,*p;

	indx=NR_ivector(1,n);
	p=NR_vector(1,n);
	fvec=NR_vector(1,n);
	fjac=NR_matrix(1,n,1,n);
	for (k=1;k<=ntrial;k++) {
		usrfun(x,n,fvec,fjac);
		errf=0.0;
		for (i=1;i<=n;i++) errf += fabs(fvec[i]);
		if (errf <= tolf) FREERETURN
		for (i=1;i<=n;i++) p[i] = -fvec[i];
		NR_ludcmp(fjac,n,indx,&d);
		NR_lubksb(fjac,n,indx,p);
		errx=0.0;
		for (i=1;i<=n;i++) {
			errx += fabs(p[i]);
			x[i] += p[i];
		}
		if (errx <= tolx) FREERETURN
	}
	FREERETURN
}
#undef FREERETURN
#undef NRANSI
