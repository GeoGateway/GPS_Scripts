#define NRANSI
#include "nr_util.h"

extern int nn;
extern float *fvec;
extern void (*nrfuncv)(int n, float v[], float f[]);

float NR_fmin(float x[])
{
	int i;
	float sum;

	(*nrfuncv)(nn,x,fvec);
	for (sum=0.0,i=1;i<=nn;i++) sum += NR_SQR(fvec[i]);
	return 0.5*sum;
}
#undef NRANSI
