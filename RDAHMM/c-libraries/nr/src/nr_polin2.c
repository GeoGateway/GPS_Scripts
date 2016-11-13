#define NRANSI
#include "nr_util.h"

void NR_polin2(float x1a[], float x2a[], float **ya, int m, int n, float x1,
	float x2, float *y, float *dy)
{
	void NR_polint(float xa[], float ya[], int n, float x, float *y, float *dy);
	int j;
	float *ymtmp;

	ymtmp=NR_vector(1,m);
	for (j=1;j<=m;j++) {
		NR_polint(x2a,ya[j],n,x2,&ymtmp[j],dy);
	}
	NR_polint(x1a,ymtmp,m,x1,y,dy);
	NR_free_vector(ymtmp,1,m);
}
#undef NRANSI
