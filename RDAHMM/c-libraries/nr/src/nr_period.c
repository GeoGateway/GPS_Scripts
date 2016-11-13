#include <math.h>
#define NRANSI
#include "nr_util.h"
#define TWOPID 6.2831853071795865

void NR_period(float x[], float y[], int n, float ofac, float hifac, float px[],
	float py[], int np, int *nout, int *jmax, float *prob)
{
	void NR_avevar(float data[], unsigned long n, float *ave, float *var);
	int i,j;
	float ave,c,cc,cwtau,effm,expy,pnow,pymax,s,ss,sumc,sumcy,sums,sumsh,
		sumsy,swtau,var,wtau,xave,xdif,xmax,xmin,yy;
	double arg,wtemp,*wi,*wpi,*wpr,*wr;

	wi=NR_dvector(1,n);
	wpi=NR_dvector(1,n);
	wpr=NR_dvector(1,n);
	wr=NR_dvector(1,n);
	*nout=0.5*ofac*hifac*n;
	if (*nout > np) NR_error("output arrays too short in NR_period");
	NR_avevar(y,n,&ave,&var);
	xmax=xmin=x[1];
	for (j=1;j<=n;j++) {
		if (x[j] > xmax) xmax=x[j];
		if (x[j] < xmin) xmin=x[j];
	}
	xdif=xmax-xmin;
	xave=0.5*(xmax+xmin);
	pymax=0.0;
	pnow=1.0/(xdif*ofac);
	for (j=1;j<=n;j++) {
		arg=TWOPID*((x[j]-xave)*pnow);
		wpr[j] = -2.0*NR_SQR(sin(0.5*arg));
		wpi[j]=sin(arg);
		wr[j]=cos(arg);
		wi[j]=wpi[j];
	}
	for (i=1;i<=(*nout);i++) {
		px[i]=pnow;
		sumsh=sumc=0.0;
		for (j=1;j<=n;j++) {
			c=wr[j];
			s=wi[j];
			sumsh += s*c;
			sumc += (c-s)*(c+s);
		}
		wtau=0.5*atan2(2.0*sumsh,sumc);
		swtau=sin(wtau);
		cwtau=cos(wtau);
		sums=sumc=sumsy=sumcy=0.0;
		for (j=1;j<=n;j++) {
			s=wi[j];
			c=wr[j];
			ss=s*cwtau-c*swtau;
			cc=c*cwtau+s*swtau;
			sums += ss*ss;
			sumc += cc*cc;
			yy=y[j]-ave;
			sumsy += yy*ss;
			sumcy += yy*cc;
			wr[j]=((wtemp=wr[j])*wpr[j]-wi[j]*wpi[j])+wr[j];
			wi[j]=(wi[j]*wpr[j]+wtemp*wpi[j])+wi[j];
		}
		py[i]=0.5*(sumcy*sumcy/sumc+sumsy*sumsy/sums)/var;
		if (py[i] >= pymax) pymax=py[(*jmax=i)];
		pnow += 1.0/(ofac*xdif);
	}
	expy=exp(-pymax);
	effm=2.0*(*nout)/ofac;
	*prob=effm*expy;
	if (*prob > 0.01) *prob=1.0-pow(1.0-expy,effm);
	NR_free_dvector(wr,1,n);
	NR_free_dvector(wpr,1,n);
	NR_free_dvector(wpi,1,n);
	NR_free_dvector(wi,1,n);
}
#undef TWOPID
#undef NRANSI
