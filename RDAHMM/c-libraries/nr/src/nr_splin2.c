#define NRANSI
#include "nr_util.h"

void NR_splin2(float x1a[], float x2a[], float **ya, float **y2a, int m, int n,
	float x1, float x2, float *y)
{
	void NR_spline(float x[], float y[], int n, float yp1, float ypn, float y2[]);
	void NR_splint(float xa[], float ya[], float y2a[], int n, float x, float *y);
	int j;
	float *ytmp,*yytmp;

	ytmp=NR_vector(1,m);
	yytmp=NR_vector(1,m);
	for (j=1;j<=m;j++)
		NR_splint(x2a,ya[j],y2a[j],n,x2,&yytmp[j]);
	NR_spline(x1a,yytmp,m,1.0e30,1.0e30,ytmp);
	NR_splint(x1a,yytmp,ytmp,m,x1,y);
	NR_free_vector(yytmp,1,m);
	NR_free_vector(ytmp,1,m);
}
#undef NRANSI
