#include <math.h>
#define NRANSI
#include "nr_util.h"
#define EPS 1.0e-10
#define FPMIN 1.0e-30
#define MAXIT 10000
#define XMIN 2.0
#ifndef PI
#ifdef M_PI
#define PI  M_PI
#else
#define PI 3.14159265358979323846
#endif
#endif

void NR_bessjy(float x, float xnu, float *NR_rj, float *ry, float *NR_rjp, float *ryp)
{
	void NR_beschb(double x, double *gam1, double *gam2, double *gampl,
		double *gammi);
	int i,isign,l,nl;
	double a,b,br,bi,c,cr,ci,d,del,del1,den,di,dlr,dli,dr,e,f,fact,fact2,
		fact3,ff,gam,gam1,gam2,gammi,gampl,h,p,pimu,pimu2,q,r,NR_rjl,
		NR_rjl1,NR_rjmu,NR_rjp1,NR_rjpl,NR_rjtemp,ry1,rymu,rymup,rytemp,sum,sum1,
		temp,w,x2,xi,xi2,xmu,xmu2;

	if (x <= 0.0 || xnu < 0.0) NR_error("bad arguments in NR_bessjy");
	nl=(x < XMIN ? (int)(xnu+0.5) : NR_IMAX(0,(int)(xnu-x+1.5)));
	xmu=xnu-nl;
	xmu2=xmu*xmu;
	xi=1.0/x;
	xi2=2.0*xi;
	w=xi2/PI;
	isign=1;
	h=xnu*xi;
	if (h < FPMIN) h=FPMIN;
	b=xi2*xnu;
	d=0.0;
	c=h;
	for (i=1;i<=MAXIT;i++) {
		b += xi2;
		d=b-d;
		if (fabs(d) < FPMIN) d=FPMIN;
		c=b-1.0/c;
		if (fabs(c) < FPMIN) c=FPMIN;
		d=1.0/d;
		del=c*d;
		h=del*h;
		if (d < 0.0) isign = -isign;
		if (fabs(del-1.0) < EPS) break;
	}
	if (i > MAXIT) NR_error("x too large in NR_bessjy; try asymptotic expansion");
	NR_rjl=isign*FPMIN;
	NR_rjpl=h*NR_rjl;
	NR_rjl1=NR_rjl;
	NR_rjp1=NR_rjpl;
	fact=xnu*xi;
	for (l=nl;l>=1;l--) {
		NR_rjtemp=fact*NR_rjl+NR_rjpl;
		fact -= xi;
		NR_rjpl=fact*NR_rjtemp-NR_rjl;
		NR_rjl=NR_rjtemp;
	}
	if (NR_rjl == 0.0) NR_rjl=EPS;
	f=NR_rjpl/NR_rjl;
	if (x < XMIN) {
		x2=0.5*x;
		pimu=PI*xmu;
		fact = (fabs(pimu) < EPS ? 1.0 : pimu/sin(pimu));
		d = -log(x2);
		e=xmu*d;
		fact2 = (fabs(e) < EPS ? 1.0 : sinh(e)/e);
		NR_beschb(xmu,&gam1,&gam2,&gampl,&gammi);
		ff=2.0/PI*fact*(gam1*cosh(e)+gam2*fact2*d);
		e=exp(e);
		p=e/(gampl*PI);
		q=1.0/(e*PI*gammi);
		pimu2=0.5*pimu;
		fact3 = (fabs(pimu2) < EPS ? 1.0 : sin(pimu2)/pimu2);
		r=PI*pimu2*fact3*fact3;
		c=1.0;
		d = -x2*x2;
		sum=ff+r*q;
		sum1=p;
		for (i=1;i<=MAXIT;i++) {
			ff=(i*ff+p+q)/(i*i-xmu2);
			c *= (d/i);
			p /= (i-xmu);
			q /= (i+xmu);
			del=c*(ff+r*q);
			sum += del;
			del1=c*p-i*del;
			sum1 += del1;
			if (fabs(del) < (1.0+fabs(sum))*EPS) break;
		}
		if (i > MAXIT) NR_error("NR_bessy series failed to converge");
		rymu = -sum;
		ry1 = -sum1*xi2;
		rymup=xmu*xi*rymu-ry1;
		NR_rjmu=w/(rymup-f*rymu);
	} else {
		a=0.25-xmu2;
		p = -0.5*xi;
		q=1.0;
		br=2.0*x;
		bi=2.0;
		fact=a*xi/(p*p+q*q);
		cr=br+q*fact;
		ci=bi+p*fact;
		den=br*br+bi*bi;
		dr=br/den;
		di = -bi/den;
		dlr=cr*dr-ci*di;
		dli=cr*di+ci*dr;
		temp=p*dlr-q*dli;
		q=p*dli+q*dlr;
		p=temp;
		for (i=2;i<=MAXIT;i++) {
			a += 2*(i-1);
			bi += 2.0;
			dr=a*dr+br;
			di=a*di+bi;
			if (fabs(dr)+fabs(di) < FPMIN) dr=FPMIN;
			fact=a/(cr*cr+ci*ci);
			cr=br+cr*fact;
			ci=bi-ci*fact;
			if (fabs(cr)+fabs(ci) < FPMIN) cr=FPMIN;
			den=dr*dr+di*di;
			dr /= den;
			di /= -den;
			dlr=cr*dr-ci*di;
			dli=cr*di+ci*dr;
			temp=p*dlr-q*dli;
			q=p*dli+q*dlr;
			p=temp;
			if (fabs(dlr-1.0)+fabs(dli) < EPS) break;
		}
		if (i > MAXIT) NR_error("cf2 failed in NR_bessjy");
		gam=(p-f)/q;
		NR_rjmu=sqrt(w/((p-f)*gam+q));
		NR_rjmu=NR_SIGN(NR_rjmu,NR_rjl);
		rymu=NR_rjmu*gam;
		rymup=rymu*(p+q/gam);
		ry1=xmu*xi*rymu-rymup;
	}
	fact=NR_rjmu/NR_rjl;
	*NR_rj=NR_rjl1*fact;
	*NR_rjp=NR_rjp1*fact;
	for (i=1;i<=nl;i++) {
		rytemp=(xmu+i)*xi2*ry1-rymu;
		rymu=ry1;
		ry1=rytemp;
	}
	*ry=rymu;
	*ryp=xnu*xi*rymu-ry1;
}
#undef EPS
#undef FPMIN
#undef MAXIT
#undef XMIN
#undef PI
#undef NRANSI
