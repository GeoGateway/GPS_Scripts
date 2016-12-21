#include <math.h>
#define NRANSI
#include "nr_util.h"

void NR_savgol(float c[], int np, int nl, int nr, int ld, int m)
{
	void NR_lubksb(float **a, int n, int *indx, float b[]);
	void NR_ludcmp(float **a, int n, int *indx, float *d);
	int imj,ipj,j,k,kk,mm,*indx;
	float d,fac,sum,**a,*b;

	if (np < nl+nr+1 || nl < 0 || nr < 0 || ld > m || nl+nr < m)
	NR_error("bad args in NR_savgol");
	indx=NR_ivector(1,m+1);
	a=NR_matrix(1,m+1,1,m+1);
	b=NR_vector(1,m+1);
	for (ipj=0;ipj<=(m << 1);ipj++) {
		sum=(ipj ? 0.0 : 1.0);
		for (k=1;k<=nr;k++) sum += pow((double)k,(double)ipj);
		for (k=1;k<=nl;k++) sum += pow((double)-k,(double)ipj);
		mm=NR_FMIN(ipj,2*m-ipj);
		for (imj = -mm;imj<=mm;imj+=2) a[1+(ipj+imj)/2][1+(ipj-imj)/2]=sum;
	}
	NR_ludcmp(a,m+1,indx,&d);
	for (j=1;j<=m+1;j++) b[j]=0.0;
	b[ld+1]=1.0;
	NR_lubksb(a,m+1,indx,b);
	for (kk=1;kk<=np;kk++) c[kk]=0.0;
	for (k = -nl;k<=nr;k++) {
		sum=b[1];
		fac=1.0;
		for (mm=1;mm<=m;mm++) sum += b[mm+1]*(fac *= k);
		kk=((np-k) % np)+1;
		c[kk]=sum;
	}
	NR_free_vector(b,1,m+1);
	NR_free_matrix(a,1,m+1,1,m+1);
	NR_free_ivector(indx,1,m+1);
}
#undef NRANSI
