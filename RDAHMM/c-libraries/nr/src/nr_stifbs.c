#include <math.h>
#define NRANSI
#include "nr_util.h"
#define KMAXX 7
#define NR_IMAXX (KMAXX+1)
#define SAFE1 0.25
#define SAFE2 0.7
#define RENR_DMAX 1.0e-5
#define RENR_DMIN 0.7
#define TINY 1.0e-30
#define SCALMX 0.1

float **d,*x;

void NR_stifbs(float y[], float dydx[], int nv, float *xx, float htry, float eps,
	float yscal[], float *hdid, float *hnext,
	void (*derivs)(float, float [], float []))
{
	void NR_jacobn(float x, float y[], float dfdx[], float **dfdy, int n);
	void NR_simpr(float y[], float dydx[], float dfdx[], float **dfdy,
		int n, float xs, float htot, int nstep, float yout[],
		void (*derivs)(float, float [], float []));
	void NR_pzextr(int iest, float xest, float yest[], float yz[], float dy[],
		int nv);
	int i,iq,k,kk,km;
	static int first=1,kmax,kopt,nvold = -1;
	static float epsold = -1.0,xnew;
	float eps1,errmax,fact,h,red,scale,work,wrkmin,xest;
	float *dfdx,**dfdy,*err,*yerr,*ysav,*yseq;
	static float a[NR_IMAXX+1];
	static float alf[KMAXX+1][KMAXX+1];
	static int nseq[NR_IMAXX+1]={0,2,6,10,14,22,34,50,70};
	int reduct,exitflag=0;

	d=NR_matrix(1,nv,1,KMAXX);
	dfdx=NR_vector(1,nv);
	dfdy=NR_matrix(1,nv,1,nv);
	err=NR_vector(1,KMAXX);
	x=NR_vector(1,KMAXX);
	yerr=NR_vector(1,nv);
	ysav=NR_vector(1,nv);
	yseq=NR_vector(1,nv);
	if(eps != epsold || nv != nvold) {
		*hnext = xnew = -1.0e29;
		eps1=SAFE1*eps;
		a[1]=nseq[1]+1;
		for (k=1;k<=KMAXX;k++) a[k+1]=a[k]+nseq[k+1];
		for (iq=2;iq<=KMAXX;iq++) {
			for (k=1;k<iq;k++)
				alf[k][iq]=pow(eps1,((a[k+1]-a[iq+1])/
					((a[iq+1]-a[1]+1.0)*(2*k+1))));
		}
		epsold=eps;
		nvold=nv;
		a[1] += nv;
		for (k=1;k<=KMAXX;k++) a[k+1]=a[k]+nseq[k+1];
		for (kopt=2;kopt<KMAXX;kopt++)
			if (a[kopt+1] > a[kopt]*alf[kopt-1][kopt]) break;
		kmax=kopt;
	}
	h=htry;
	for (i=1;i<=nv;i++) ysav[i]=y[i];
	NR_jacobn(*xx,y,dfdx,dfdy,nv);
	if (*xx != xnew || h != (*hnext)) {
		first=1;
		kopt=kmax;
	}
	reduct=0;
	for (;;) {
		for (k=1;k<=kmax;k++) {
			xnew=(*xx)+h;
			if (xnew == (*xx)) NR_error("step size underflow in NR_stifbs");
			NR_simpr(ysav,dydx,dfdx,dfdy,nv,*xx,h,nseq[k],yseq,derivs);
			xest=NR_SQR(h/nseq[k]);
			NR_pzextr(k,xest,yseq,y,yerr,nv);
			if (k != 1) {
				errmax=TINY;
				for (i=1;i<=nv;i++) errmax=NR_FMAX(errmax,fabs(yerr[i]/yscal[i]));
				errmax /= eps;
				km=k-1;
				err[km]=pow(errmax/SAFE1,1.0/(2*km+1));
			}
			if (k != 1 && (k >= kopt-1 || first)) {
				if (errmax < 1.0) {
					exitflag=1;
					break;
				}
				if (k == kmax || k == kopt+1) {
					red=SAFE2/err[km];
					break;
				}
				else if (k == kopt && alf[kopt-1][kopt] < err[km]) {
						red=1.0/err[km];
						break;
					}
				else if (kopt == kmax && alf[km][kmax-1] < err[km]) {
						red=alf[km][kmax-1]*SAFE2/err[km];
						break;
					}
				else if (alf[km][kopt] < err[km]) {
					red=alf[km][kopt-1]/err[km];
					break;
				}
			}
		}
		if (exitflag) break;
		red=NR_FMIN(red,RENR_DMIN);
		red=NR_FMAX(red,RENR_DMAX);
		h *= red;
		reduct=1;
	}
	*xx=xnew;
	*hdid=h;
	first=0;
	wrkmin=1.0e35;
	for (kk=1;kk<=km;kk++) {
		fact=NR_FMAX(err[kk],SCALMX);
		work=fact*a[kk+1];
		if (work < wrkmin) {
			scale=fact;
			wrkmin=work;
			kopt=kk+1;
		}
	}
	*hnext=h/scale;
	if (kopt >= k && kopt != kmax && !reduct) {
		fact=NR_FMAX(scale/alf[kopt-1][kopt],SCALMX);
		if (a[kopt+1]*fact <= wrkmin) {
			*hnext=h/fact;
			kopt++;
		}
	}
	NR_free_vector(yseq,1,nv);
	NR_free_vector(ysav,1,nv);
	NR_free_vector(yerr,1,nv);
	NR_free_vector(x,1,KMAXX);
	NR_free_vector(err,1,KMAXX);
	NR_free_matrix(dfdy,1,nv,1,nv);
	NR_free_vector(dfdx,1,nv);
	NR_free_matrix(d,1,nv,1,KMAXX);
}
#undef KMAXX
#undef NR_IMAXX
#undef SAFE1
#undef SAFE2
#undef RENR_DMAX
#undef RENR_DMIN
#undef TINY
#undef SCALMX
#undef NRANSI
