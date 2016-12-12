#define NRANSI
#include "nr_util.h"

void NR_fred2(int n, float a, float b, float t[], float f[], float w[],
	float (*g)(float), float (*ak)(float, float))
{
	void NR_gauleg(float x1, float x2, float x[], float w[], int n);
	void NR_lubksb(float **a, int n, int *indx, float b[]);
	void NR_ludcmp(float **a, int n, int *indx, float *d);
	int i,j,*indx;
	float d,**omk;

	indx=NR_ivector(1,n);
	omk=NR_matrix(1,n,1,n);
	NR_gauleg(a,b,t,w,n);
	for (i=1;i<=n;i++) {
		for (j=1;j<=n;j++)
			omk[i][j]=(float)(i == j)-(*ak)(t[i],t[j])*w[j];
		f[i]=(*g)(t[i]);
	}
	NR_ludcmp(omk,n,indx,&d);
	NR_lubksb(omk,n,indx,f);
	NR_free_matrix(omk,1,n,1,n);
	NR_free_ivector(indx,1,n);
}
#undef NRANSI
