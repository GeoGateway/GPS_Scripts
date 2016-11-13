#define NRANSI
#include "nr_util.h"
#define TOL 2.0e-4

int ncom;
float *pcom,*xicom,(*nrfunc)(float []);

void NR_linmin(float p[], float xi[], int n, float *fret, float (*func)(float []))
{
	float NR_brent(float ax, float bx, float cx,
		float (*f)(float), float tol, float *xmin);
	float NR_f1dim(float x);
	void NR_mnbrak(float *ax, float *bx, float *cx, float *fa, float *fb,
		float *fc, float (*func)(float));
	int j;
	float xx,xmin,fx,fb,fa,bx,ax;

	ncom=n;
	pcom=NR_vector(1,n);
	xicom=NR_vector(1,n);
	nrfunc=func;
	for (j=1;j<=n;j++) {
		pcom[j]=p[j];
		xicom[j]=xi[j];
	}
	ax=0.0;
	xx=1.0;
	NR_mnbrak(&ax,&xx,&bx,&fa,&fx,&fb,NR_f1dim);
	*fret=NR_brent(ax,xx,bx,NR_f1dim,TOL,&xmin);
	for (j=1;j<=n;j++) {
		xi[j] *= xmin;
		p[j] += xi[j];
	}
	NR_free_vector(xicom,1,n);
	NR_free_vector(pcom,1,n);
}
#undef TOL
#undef NRANSI
