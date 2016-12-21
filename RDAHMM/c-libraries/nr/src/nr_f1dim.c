#define NRANSI
#include "nr_util.h"

extern int ncom;
extern float *pcom,*xicom,(*nrfunc)(float []);

float NR_f1dim(float x)
{
	int j;
	float f,*xt;

	xt=NR_vector(1,ncom);
	for (j=1;j<=ncom;j++) xt[j]=pcom[j]+x*xicom[j];
	f=(*nrfunc)(xt);
	NR_free_vector(xt,1,ncom);
	return f;
}
#undef NRANSI
