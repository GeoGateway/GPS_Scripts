#define NRANSI
#include "nr_util.h"

extern int ncom;
extern float *pcom,*xicom,(*nrfunc)(float []);
extern void (*nrdfun)(float [], float []);

float NR_df1dim(float x)
{
	int j;
	float df1=0.0;
	float *xt,*df;

	xt=NR_vector(1,ncom);
	df=NR_vector(1,ncom);
	for (j=1;j<=ncom;j++) xt[j]=pcom[j]+x*xicom[j];
	(*nrdfun)(xt,df);
	for (j=1;j<=ncom;j++) df1 += df[j]*xicom[j];
	NR_free_vector(df,1,ncom);
	NR_free_vector(xt,1,ncom);
	return df1;
}
#undef NRANSI
