#include <math.h>
#define NRANSI
#include "nr_util.h"

void NR_ksone(float data[], unsigned long n, float (*func)(float), float *d,
	float *prob)
{
	float NR_probks(float alam);
	void NR_sort(unsigned long n, float arr[]);
	unsigned long j;
	float dt,en,ff,fn,fo=0.0;

	NR_sort(n,data);
	en=n;
	*d=0.0;
	for (j=1;j<=n;j++) {
		fn=j/en;
		ff=(*func)(data[j]);
		dt=NR_FMAX(fabs(fo-ff),fabs(fn-ff));
		if (dt > *d) *d=dt;
		fo=fn;
	}
	en=sqrt(en);
	*prob=NR_probks((en+0.12+0.11/en)*(*d));
}
#undef NRANSI
