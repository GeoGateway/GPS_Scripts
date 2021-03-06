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

double x;

void NR_quadmx(float **a, int n)
{
	void NR_kermom(double w[], double y, int m);
	void NR_wwghts(float wghts[], int n, float h,
		void (*NR_kermom)(double [], double ,int));
	int j,k;
	float h,*wt,xx,cx;

	wt=NR_vector(1,n);
	h=PI/(n-1);
	for (j=1;j<=n;j++) {
		x=xx=(j-1)*h;
		NR_wwghts(wt,n,h,NR_kermom);
		cx=cos(xx);
		for (k=1;k<=n;k++) a[j][k]=wt[k]*cx*cos((k-1)*h);
		++a[j][j];
	}
	NR_free_vector(wt,1,n);
}
#undef PI
#undef NRANSI
