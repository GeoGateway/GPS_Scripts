#define NRANSI
#include "nr_util.h"

void NR_mrqmin(float x[], float y[], float sig[], int ndata, float a[], int ia[],
	int ma, float **covar, float **alpha, float *chisq,
	void (*funcs)(float, float [], float *, float [], int), float *alamda)
{
	void NR_covsrt(float **covar, int ma, int ia[], int mNR_fit);
	void NR_gaussj(float **a, int n, float **b, int m);
	void NR_mrqcof(float x[], float y[], float sig[], int ndata, float a[],
		int ia[], int ma, float **alpha, float NR_beta[], float *chisq,
		void (*funcs)(float, float [], float *, float [], int));
	int j,k,l;
	static int mNR_fit;
	static float ochisq,*atry,*NR_beta,*da,**oneda;

	if (*alamda < 0.0) {
		atry=NR_vector(1,ma);
		NR_beta=NR_vector(1,ma);
		da=NR_vector(1,ma);
		for (mNR_fit=0,j=1;j<=ma;j++)
			if (ia[j]) mNR_fit++;
		oneda=NR_matrix(1,mNR_fit,1,1);
		*alamda=0.001;
		NR_mrqcof(x,y,sig,ndata,a,ia,ma,alpha,NR_beta,chisq,funcs);
		ochisq=(*chisq);
		for (j=1;j<=ma;j++) atry[j]=a[j];
	}
	for (j=1;j<=mNR_fit;j++) {
		for (k=1;k<=mNR_fit;k++) covar[j][k]=alpha[j][k];
		covar[j][j]=alpha[j][j]*(1.0+(*alamda));
		oneda[j][1]=NR_beta[j];
	}
	NR_gaussj(covar,mNR_fit,oneda,1);
	for (j=1;j<=mNR_fit;j++) da[j]=oneda[j][1];
	if (*alamda == 0.0) {
		NR_covsrt(covar,ma,ia,mNR_fit);
		NR_free_matrix(oneda,1,mNR_fit,1,1);
		NR_free_vector(da,1,ma);
		NR_free_vector(NR_beta,1,ma);
		NR_free_vector(atry,1,ma);
		return;
	}
	for (j=0,l=1;l<=ma;l++)
		if (ia[l]) atry[l]=a[l]+da[++j];
	NR_mrqcof(x,y,sig,ndata,atry,ia,ma,covar,da,chisq,funcs);
	if (*chisq < ochisq) {
		*alamda *= 0.1;
		ochisq=(*chisq);
		for (j=1;j<=mNR_fit;j++) {
			for (k=1;k<=mNR_fit;k++) alpha[j][k]=covar[j][k];
			NR_beta[j]=da[j];
		}
		for (l=1;l<=ma;l++) a[l]=atry[l];
	} else {
		*alamda *= 10.0;
		*chisq=ochisq;
	}
}
#undef NRANSI
