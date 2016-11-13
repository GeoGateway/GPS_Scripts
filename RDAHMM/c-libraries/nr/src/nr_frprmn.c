#include <math.h>
#define NRANSI
#include "nr_util.h"
#define ITMAX 200
#define EPS 1.0e-10
#define FREEALL NR_free_vector(xi,1,n);NR_free_vector(h,1,n);NR_free_vector(g,1,n);

void NR_frprmn(float p[], int n, float ftol, int *iter, float *fret,
	float (*func)(float []), void (*dfunc)(float [], float []))
{
	void NR_linmin(float p[], float xi[], int n, float *fret,
		float (*func)(float []));
	int j,its;
	float gg,gam,fp,dgg;
	float *g,*h,*xi;

	g=NR_vector(1,n);
	h=NR_vector(1,n);
	xi=NR_vector(1,n);
	fp=(*func)(p);
	(*dfunc)(p,xi);
	for (j=1;j<=n;j++) {
		g[j] = -xi[j];
		xi[j]=h[j]=g[j];
	}
	for (its=1;its<=ITMAX;its++) {
		*iter=its;
		NR_linmin(p,xi,n,fret,func);
		if (2.0*fabs(*fret-fp) <= ftol*(fabs(*fret)+fabs(fp)+EPS)) {
			FREEALL
			return;
		}
		fp=(*func)(p);
		(*dfunc)(p,xi);
		dgg=gg=0.0;
		for (j=1;j<=n;j++) {
			gg += g[j]*g[j];
			dgg += (xi[j]+g[j])*xi[j];
		}
		if (gg == 0.0) {
			FREEALL
			return;
		}
		gam=dgg/gg;
		for (j=1;j<=n;j++) {
			g[j] = -xi[j];
			xi[j]=h[j]=g[j]+gam*h[j];
		}
	}
	NR_error("Too many iterations in NR_frprmn");
}
#undef ITMAX
#undef EPS
#undef FREEALL
#undef NRANSI
