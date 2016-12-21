#define NRANSI
#include "nr_util.h"

void lNR_fit(float x[], float y[], float sig[], int ndat, float a[], int ia[],
	int ma, float **covar, float *chisq, void (*funcs)(float, float [], int))
{
	void NR_covsrt(float **covar, int ma, int ia[], int mNR_fit);
	void NR_gaussj(float **a, int n, float **b, int m);
	int i,j,k,l,m,mNR_fit=0;
	float ym,wt,sum,sig2i,**NR_beta,*afunc;

	NR_beta=NR_matrix(1,ma,1,1);
	afunc=NR_vector(1,ma);
	for (j=1;j<=ma;j++)
		if (ia[j]) mNR_fit++;
	if (mNR_fit == 0) NR_error("lNR_fit: no parameters to be NR_fitted");
	for (j=1;j<=mNR_fit;j++) {
		for (k=1;k<=mNR_fit;k++) covar[j][k]=0.0;
		NR_beta[j][1]=0.0;
	}
	for (i=1;i<=ndat;i++) {
		(*funcs)(x[i],afunc,ma);
		ym=y[i];
		if (mNR_fit < ma) {
			for (j=1;j<=ma;j++)
				if (!ia[j]) ym -= a[j]*afunc[j];
		}
		sig2i=1.0/NR_SQR(sig[i]);
		for (j=0,l=1;l<=ma;l++) {
			if (ia[l]) {
				wt=afunc[l]*sig2i;
				for (j++,k=0,m=1;m<=l;m++)
					if (ia[m]) covar[j][++k] += wt*afunc[m];
				NR_beta[j][1] += ym*wt;
			}
		}
	}
	for (j=2;j<=mNR_fit;j++)
		for (k=1;k<j;k++)
			covar[k][j]=covar[j][k];
	NR_gaussj(covar,mNR_fit,NR_beta,1);
	for (j=0,l=1;l<=ma;l++)
		if (ia[l]) a[l]=NR_beta[++j][1];
	*chisq=0.0;
	for (i=1;i<=ndat;i++) {
		(*funcs)(x[i],afunc,ma);
		for (sum=0.0,j=1;j<=ma;j++) sum += a[j]*afunc[j];
		*chisq += NR_SQR((y[i]-sum)/sig[i]);
	}
	NR_covsrt(covar,ma,ia,mNR_fit);
	NR_free_vector(afunc,1,ma);
	NR_free_matrix(NR_beta,1,ma,1,1);
}
#undef NRANSI
