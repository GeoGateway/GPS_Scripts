#define NRANSI
#include "nr_util.h"
#define TOL 1.0e-5

void NR_svdfit(float x[], float y[], float sig[], int ndata, float a[], int ma,
	float **u, float **v, float w[], float *chisq,
	void (*funcs)(float, float [], int))
{
	void NR_svbksb(float **u, float w[], float **v, int m, int n, float b[],
		float x[]);
	void NR_svdcmp(float **a, int m, int n, float w[], float **v);
	int j,i;
	float wmax,tmp,thresh,sum,*b,*afunc;

	b=NR_vector(1,ndata);
	afunc=NR_vector(1,ma);
	for (i=1;i<=ndata;i++) {
		(*funcs)(x[i],afunc,ma);
		tmp=1.0/sig[i];
		for (j=1;j<=ma;j++) u[i][j]=afunc[j]*tmp;
		b[i]=y[i]*tmp;
	}
	NR_svdcmp(u,ndata,ma,w,v);
	wmax=0.0;
	for (j=1;j<=ma;j++)
		if (w[j] > wmax) wmax=w[j];
	thresh=TOL*wmax;
	for (j=1;j<=ma;j++)
		if (w[j] < thresh) w[j]=0.0;
	NR_svbksb(u,w,v,ndata,ma,b,a);
	*chisq=0.0;
	for (i=1;i<=ndata;i++) {
		(*funcs)(x[i],afunc,ma);
		for (sum=0.0,j=1;j<=ma;j++) sum += a[j]*afunc[j];
		*chisq += (tmp=(y[i]-sum)/sig[i],tmp*tmp);
	}
	NR_free_vector(afunc,1,ma);
	NR_free_vector(b,1,ndata);
}
#undef TOL
#undef NRANSI
