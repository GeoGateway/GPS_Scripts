#include <math.h>
#define NRANSI
#include "nr_util.h"
#define EPS 1.0e-7

extern int ndatat;
extern float *xt,*yt,aa,abdevt;

float NR_rofunc(float b)
{
	float NR_select(unsigned long k, unsigned long n, float arr[]);
	int j;
	float *arr,d,sum=0.0;

	arr=NR_vector(1,ndatat);
	for (j=1;j<=ndatat;j++) arr[j]=yt[j]-b*xt[j];
	if (ndatat & 1) {
		aa=NR_select((ndatat+1)>>1,ndatat,arr);
	}
	else {
		j=ndatat >> 1;
		aa=0.5*(NR_select(j,ndatat,arr)+NR_select(j+1,ndatat,arr));
	}
	abdevt=0.0;
	for (j=1;j<=ndatat;j++) {
		d=yt[j]-(b*xt[j]+aa);
		abdevt += fabs(d);
		if (yt[j] != 0.0) d /= fabs(yt[j]);
		if (fabs(d) > EPS) sum += (d >= 0.0 ? xt[j] : -xt[j]);
	}
	NR_free_vector(arr,1,ndatat);
	return sum;
}
#undef EPS
#undef NRANSI
