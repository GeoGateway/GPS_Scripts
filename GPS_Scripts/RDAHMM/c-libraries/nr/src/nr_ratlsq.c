#include <stdio.h>
#include <math.h>
#define NRANSI
#include "nr_util.h"
#define NPFAC 8
#define MAXIT 5
#define PIO2 (3.141592653589793/2.0)
#define BIG 1.0e30

void NR_ratlsq(double (*fn)(double), double a, double b, int mm, int kk,
	double cof[], double *dev)
{
	double NR_ratval(double x, double cof[], int mm, int kk);
	void NR_dsvbksb(double **u, double w[], double **v, int m, int n, double b[],
		double x[]);
	void NR_dsvdcmp(double **a, int m, int n, double w[], double **v);
	int i,it,j,ncof,npt;
	double devmax,e,hth,power,sum,*bb,*coff,*ee,*fs,**u,**v,*w,*wt,*xs;

	ncof=mm+kk+1;
	npt=NPFAC*ncof;
	bb=NR_dvector(1,npt);
	coff=NR_dvector(0,ncof-1);
	ee=NR_dvector(1,npt);
	fs=NR_dvector(1,npt);
	u=NR_dmatrix(1,npt,1,ncof);
	v=NR_dmatrix(1,ncof,1,ncof);
	w=NR_dvector(1,ncof);
	wt=NR_dvector(1,npt);
	xs=NR_dvector(1,npt);
	*dev=BIG;
	for (i=1;i<=npt;i++) {
		if (i < npt/2) {
			hth=PIO2*(i-1)/(npt-1.0);
			xs[i]=a+(b-a)*NR_DSQR(sin(hth));
		} else {
			hth=PIO2*(npt-i)/(npt-1.0);
			xs[i]=b-(b-a)*NR_DSQR(sin(hth));
		}
		fs[i]=(*fn)(xs[i]);
		wt[i]=1.0;
		ee[i]=1.0;
	}
	e=0.0;
	for (it=1;it<=MAXIT;it++) {
		for (i=1;i<=npt;i++) {
			power=wt[i];
			bb[i]=power*(fs[i]+NR_SIGN(e,ee[i]));
			for (j=1;j<=mm+1;j++) {
				u[i][j]=power;
				power *= xs[i];
			}
			power = -bb[i];
			for (j=mm+2;j<=ncof;j++) {
				power *= xs[i];
				u[i][j]=power;
			}
		}
		NR_dsvdcmp(u,npt,ncof,w,v);
		NR_dsvbksb(u,w,v,npt,ncof,bb,coff-1);
		devmax=sum=0.0;
		for (j=1;j<=npt;j++) {
			ee[j]=NR_ratval(xs[j],coff,mm,kk)-fs[j];
			wt[j]=fabs(ee[j]);
			sum += wt[j];
			if (wt[j] > devmax) devmax=wt[j];
		}
		e=sum/npt;
		if (devmax <= *dev) {
			for (j=0;j<ncof;j++) cof[j]=coff[j];
			*dev=devmax;
		}
		printf(" NR_ratlsq iteration= %2d  max error= %10.3e\n",it,devmax);
	}
	NR_free_dvector(xs,1,npt);
	NR_free_dvector(wt,1,npt);
	NR_free_dvector(w,1,ncof);
	NR_free_dmatrix(v,1,ncof,1,ncof);
	NR_free_dmatrix(u,1,npt,1,ncof);
	NR_free_dvector(fs,1,npt);
	NR_free_dvector(ee,1,npt);
	NR_free_dvector(coff,0,ncof-1);
	NR_free_dvector(bb,1,npt);
}
#undef NPFAC
#undef MAXIT
#undef PIO2
#undef BIG
#undef NRANSI
