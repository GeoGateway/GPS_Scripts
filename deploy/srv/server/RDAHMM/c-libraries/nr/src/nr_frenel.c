#include <math.h>
#include "nr_complex.h"
#define EPS 6.0e-8
#define MAXIT 100
#define FPMIN 1.0e-30
#define XMIN 1.5
#ifndef PI
#ifdef M_PI
#define PI  M_PI
#else
#define PI 3.14159265358979323846
#endif
#endif
#define PIBY2 (PI/2.0)
#define TRUE 1
#define ONE NR_Complex(1.0,0.0)

void NR_frenel(float x, float *s, float *c)
{
	void NR_error(char error_text[]);
	int k,n,odd;
	float a,ax,fact,pix2,sign,sum,sumc,sums,term,test;
	nr_fcomplex b,cc,d,h,del,cs;

	ax=fabs(x);
	if (ax < sqrt(FPMIN)) {
		*s=0.0;
		*c=ax;
	} else if (ax <= XMIN) {
		sum=sums=0.0;
		sumc=ax;
		sign=1.0;
		fact=PIBY2*ax*ax;
		odd=TRUE;
		term=ax;
		n=3;
		for (k=1;k<=MAXIT;k++) {
			term *= fact/k;
			sum += sign*term/n;
			test=fabs(sum)*EPS;
			if (odd) {
				sign = -sign;
				sums=sum;
				sum=sumc;
			} else {
				sumc=sum;
				sum=sums;
			}
			if (term < test) break;
			odd=!odd;
			n += 2;
		}
		if (k > MAXIT) NR_error("series failed in NR_frenel");
		*s=sums;
		*c=sumc;
	} else {
		pix2=PI*ax*ax;
		b=NR_Complex(1.0,-pix2);
		cc=NR_Complex(1.0/FPMIN,0.0);
		d=h=NR_Cdiv(ONE,b);
		n = -1;
		for (k=2;k<=MAXIT;k++) {
			n += 2;
			a = -n*(n+1);
			b=NR_Cadd(b,NR_Complex(4.0,0.0));
			d=NR_Cdiv(ONE,NR_Cadd(NR_RCmul(a,d),b));
			cc=NR_Cadd(b,NR_Cdiv(NR_Complex(a,0.0),cc));
			del=NR_Cmul(cc,d);
			h=NR_Cmul(h,del);
			if (fabs(del.r-1.0)+fabs(del.i) < EPS) break;
		}
		if (k > MAXIT) NR_error("cf failed in NR_frenel");
		h=NR_Cmul(NR_Complex(ax,-ax),h);
		cs=NR_Cmul(NR_Complex(0.5,0.5),
			NR_Csub(ONE,NR_Cmul(NR_Complex(cos(0.5*pix2),sin(0.5*pix2)),h)));
		*c=cs.r;
		*s=cs.i;
	}
	if (x < 0.0) {
		*c = -(*c);
		*s = -(*s);
	}
}
#undef EPS
#undef MAXIT
#undef FPMIN
#undef XMIN
#undef PI
#undef PIBY2
#undef TRUE
#undef ONE
