#include <math.h>
#define NRANSI
#include "nr_util.h"
#define ITMAX 200

void NR_powell(float p[], float **xi, int n, float ftol, int *iter, float *fret,
	float (*func)(float []))
{
	void NR_linmin(float p[], float xi[], int n, float *fret,
		float (*func)(float []));
	int i,ibig,j;
	float del,fp,fptt,t,*pt,*ptt,*xit;

	pt=NR_vector(1,n);
	ptt=NR_vector(1,n);
	xit=NR_vector(1,n);
	*fret=(*func)(p);
	for (j=1;j<=n;j++) pt[j]=p[j];
	for (*iter=1;;++(*iter)) {
		fp=(*fret);
		ibig=0;
		del=0.0;
		for (i=1;i<=n;i++) {
			for (j=1;j<=n;j++) xit[j]=xi[j][i];
			fptt=(*fret);
			NR_linmin(p,xit,n,fret,func);
			if (fabs(fptt-(*fret)) > del) {
				del=fabs(fptt-(*fret));
				ibig=i;
			}
		}
		if (2.0*fabs(fp-(*fret)) <= ftol*(fabs(fp)+fabs(*fret))) {
			NR_free_vector(xit,1,n);
			NR_free_vector(ptt,1,n);
			NR_free_vector(pt,1,n);
			return;
		}
		if (*iter == ITMAX) NR_error("NR_powell exceeding maximum iterations.");
		for (j=1;j<=n;j++) {
			ptt[j]=2.0*p[j]-pt[j];
			xit[j]=p[j]-pt[j];
			pt[j]=p[j];
		}
		fptt=(*func)(ptt);
		if (fptt < fp) {
			t=2.0*(fp-2.0*(*fret)+fptt)*NR_SQR(fp-(*fret)-del)-del*NR_SQR(fp-fptt);
			if (t < 0.0) {
				NR_linmin(p,xit,n,fret,func);
				for (j=1;j<=n;j++) {
					xi[j][ibig]=xi[j][n];
					xi[j][n]=xit[j];
				}
			}
		}
	}
}
#undef ITMAX
#undef NRANSI
