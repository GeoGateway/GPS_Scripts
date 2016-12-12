#define NRANSI
#include "nr_util.h"

float **y,*xx;

void NR_rkdumb(float vstart[], int nvar, float x1, float x2, int nstep,
	void (*derivs)(float, float [], float []))
{
	void NR_rk4(float y[], float dydx[], int n, float x, float h, float yout[],
		void (*derivs)(float, float [], float []));
	int i,k;
	float x,h;
	float *v,*vout,*dv;

	v=NR_vector(1,nvar);
	vout=NR_vector(1,nvar);
	dv=NR_vector(1,nvar);
	for (i=1;i<=nvar;i++) {
		v[i]=vstart[i];
		y[i][1]=v[i];
	}
	xx[1]=x1;
	x=x1;
	h=(x2-x1)/nstep;
	for (k=1;k<=nstep;k++) {
		(*derivs)(x,v,dv);
		NR_rk4(v,dv,nvar,x,h,vout,derivs);
		if ((float)(x+h) == x) NR_error("Step size too small in routine NR_rkdumb");
		x += h;
		xx[k+1]=x;
		for (i=1;i<=nvar;i++) {
			v[i]=vout[i];
			y[i][k+1]=v[i];
		}
	}
	NR_free_vector(dv,1,nvar);
	NR_free_vector(vout,1,nvar);
	NR_free_vector(v,1,nvar);
}
#undef NRANSI
