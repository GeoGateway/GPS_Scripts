#define NRANSI
#include "nr_util.h"

void NR_mrqcof(float x[], float y[], float sig[], int ndata, float a[], int ia[],
	int ma, float **alpha, float NR_beta[], float *chisq,
	void (*funcs)(float, float [], float *, float [], int))
{
	int i,j,k,l,m,mNR_fit=0;
	float ymod,wt,sig2i,dy,*dyda;

	dyda=NR_vector(1,ma);
	for (j=1;j<=ma;j++)
		if (ia[j]) mNR_fit++;
	for (j=1;j<=mNR_fit;j++) {
		for (k=1;k<=j;k++) alpha[j][k]=0.0;
		NR_beta[j]=0.0;
	}
	*chisq=0.0;
	for (i=1;i<=ndata;i++) {
		(*funcs)(x[i],a,&ymod,dyda,ma);
		sig2i=1.0/(sig[i]*sig[i]);
		dy=y[i]-ymod;
		for (j=0,l=1;l<=ma;l++) {
			if (ia[l]) {
				wt=dyda[l]*sig2i;
				for (j++,k=0,m=1;m<=l;m++)
					if (ia[m]) alpha[j][++k] += wt*dyda[m];
				NR_beta[j] += dy*wt;
			}
		}
		*chisq += dy*dy*sig2i;
	}
	for (j=2;j<=mNR_fit;j++)
		for (k=1;k<j;k++) alpha[k][j]=alpha[j][k];
	NR_free_vector(dyda,1,ma);
}
#undef NRANSI
