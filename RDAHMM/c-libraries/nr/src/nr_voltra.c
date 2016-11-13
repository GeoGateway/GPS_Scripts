#define NRANSI
#include "nr_util.h"

void NR_voltra(int n, int m, float t0, float h, float *t, float **f,
	float (*g)(int, float), float (*ak)(int, int, float, float))
{
	void NR_lubksb(float **a, int n, int *indx, float b[]);
	void NR_ludcmp(float **a, int n, int *indx, float *d);
	int i,j,k,l,*indx;
	float d,sum,**a,*b;

	indx=NR_ivector(1,m);
	a=NR_matrix(1,m,1,m);
	b=NR_vector(1,m);
	t[1]=t0;
	for (k=1;k<=m;k++) f[k][1]=(*g)(k,t[1]);
	for (i=2;i<=n;i++) {
		t[i]=t[i-1]+h;
		for (k=1;k<=m;k++) {
			sum=(*g)(k,t[i]);
			for (l=1;l<=m;l++) {
				sum += 0.5*h*(*ak)(k,l,t[i],t[1])*f[l][1];
				for (j=2;j<i;j++)
					sum += h*(*ak)(k,l,t[i],t[j])*f[l][j];
				a[k][l]=(k == l)-0.5*h*(*ak)(k,l,t[i],t[i]);
			}
			b[k]=sum;
		}
		NR_ludcmp(a,m,indx,&d);
		NR_lubksb(a,m,indx,b);
		for (k=1;k<=m;k++) f[k][i]=b[k];
	}
	NR_free_vector(b,1,m);
	NR_free_matrix(a,1,m,1,m);
	NR_free_ivector(indx,1,m);
}
#undef NRANSI
