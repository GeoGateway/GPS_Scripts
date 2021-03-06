#include <stdio.h>
#include <math.h>
#define NRANSI
#include "nr_util.h"
#ifndef PI
#ifdef M_PI
#define PI  M_PI
#else
#define PI 3.14159265358979323846
#endif
#endif
#define N 40

int main(void)  /* Program fredex */
{
	void NR_lubksb(float **a, int n, int *indx, float b[]);
	void NR_ludcmp(float **a, int n, int *indx, float *d);
	void NR_quadmx(float **a, int n);
	float **a,d,*g,x;
	int *indx,j;

	indx=NR_ivector(1,N);
	a=NR_matrix(1,N,1,N);
	g=NR_vector(1,N);
	NR_quadmx(a,N);
	NR_ludcmp(a,N,indx,&d);
	for (j=1;j<=N;j++) g[j]=sin((j-1)*PI/(N-1));
	NR_lubksb(a,N,indx,g);
	for (j=1;j<=N;j++) {
		x=(j-1)*PI/(N-1);
		printf("%6.2d %12.6f %12.6f\n",j,x,g[j]);
	}
	NR_free_vector(g,1,N);
	NR_free_matrix(a,1,N,1,N);
	NR_free_ivector(indx,1,N);
	return 0;
}
#undef N
#undef PI
#undef NRANSI
