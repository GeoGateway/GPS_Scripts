#include <math.h>
#define NRANSI
#include "nr_util.h"
#define MAXSTP 10000
#define TINY 1.0e-30

extern int kmax,kount;
extern float *xp,**yp,dxsav;

void NR_odeint(float ystart[], int nvar, float x1, float x2, float eps, float h1,
	float hmin, int *nok, int *nbad,
	void (*derivs)(float, float [], float []),
	void (*NR_rkqs)(float [], float [], int, float *, float, float, float [],
	float *, float *, void (*)(float, float [], float [])))
{
	int nstp,i;
	float xsav,x,hnext,hdid,h;
	float *yscal,*y,*dydx;

	yscal=NR_vector(1,nvar);
	y=NR_vector(1,nvar);
	dydx=NR_vector(1,nvar);
	x=x1;
	h=NR_SIGN(h1,x2-x1);
	*nok = (*nbad) = kount = 0;
	for (i=1;i<=nvar;i++) y[i]=ystart[i];
	if (kmax > 0) xsav=x-dxsav*2.0;
	for (nstp=1;nstp<=MAXSTP;nstp++) {
		(*derivs)(x,y,dydx);
		for (i=1;i<=nvar;i++)
			yscal[i]=fabs(y[i])+fabs(dydx[i]*h)+TINY;
		if (kmax > 0 && kount < kmax-1 && fabs(x-xsav) > fabs(dxsav)) {
			xp[++kount]=x;
			for (i=1;i<=nvar;i++) yp[i][kount]=y[i];
			xsav=x;
		}
		if ((x+h-x2)*(x+h-x1) > 0.0) h=x2-x;
		(*NR_rkqs)(y,dydx,nvar,&x,h,eps,yscal,&hdid,&hnext,derivs);
		if (hdid == h) ++(*nok); else ++(*nbad);
		if ((x-x2)*(x2-x1) >= 0.0) {
			for (i=1;i<=nvar;i++) ystart[i]=y[i];
			if (kmax) {
				xp[++kount]=x;
				for (i=1;i<=nvar;i++) yp[i][kount]=y[i];
			}
			NR_free_vector(dydx,1,nvar);
			NR_free_vector(y,1,nvar);
			NR_free_vector(yscal,1,nvar);
			return;
		}
		if (fabs(hnext) <= hmin) NR_error("Step size too small in NR_odeint");
		h=hnext;
	}
	NR_error("Too many steps in routine NR_odeint");
}
#undef MAXSTP
#undef TINY
#undef NRANSI
