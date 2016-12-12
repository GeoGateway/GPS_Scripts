#include <math.h>
#define NRANSI
#include "nr_util.h"
#define SAFETY 0.9
#define PGROW -0.2
#define PSHRNK -0.25
#define ERRCON 1.89e-4

void NR_rkqs(float y[], float dydx[], int n, float *x, float htry, float eps,
	float yscal[], float *hdid, float *hnext,
	void (*derivs)(float, float [], float []))
{
	void NR_rkck(float y[], float dydx[], int n, float x, float h,
		float yout[], float yerr[], void (*derivs)(float, float [], float []));
	int i;
	float errmax,h,htemp,xnew,*yerr,*ytemp;

	yerr=NR_vector(1,n);
	ytemp=NR_vector(1,n);
	h=htry;
	for (;;) {
		NR_rkck(y,dydx,n,*x,h,ytemp,yerr,derivs);
		errmax=0.0;
		for (i=1;i<=n;i++) errmax=NR_FMAX(errmax,fabs(yerr[i]/yscal[i]));
		errmax /= eps;
		if (errmax <= 1.0) break;
		htemp=SAFETY*h*pow(errmax,PSHRNK);
		h=(h >= 0.0 ? NR_FMAX(htemp,0.1*h) : NR_FMIN(htemp,0.1*h));
		xnew=(*x)+h;
		if (xnew == *x) NR_error("stepsize underflow in NR_rkqs");
	}
	if (errmax > ERRCON) *hnext=SAFETY*h*pow(errmax,PGROW);
	else *hnext=5.0*h;
	*x += (*hdid=h);
	for (i=1;i<=n;i++) y[i]=ytemp[i];
	NR_free_vector(ytemp,1,n);
	NR_free_vector(yerr,1,n);
}
#undef SAFETY
#undef PGROW
#undef PSHRNK
#undef ERRCON
#undef NRANSI
