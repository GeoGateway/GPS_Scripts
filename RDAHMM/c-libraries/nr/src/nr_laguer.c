#include <math.h>
#include "nr_complex.h"
#define NRANSI
#include "nr_util.h"
#define EPSS 1.0e-7
#define MR 8
#define MT 10
#define MAXIT (MT*MR)

void NR_laguer(nr_fcomplex a[], int m, nr_fcomplex *x, int *its)
{
	int iter,j;
	float abx,abp,abm,err;
	nr_fcomplex dx,x1,b,d,f,g,h,sq,gp,gm,g2;
	static float frac[MR+1] = {0.0,0.5,0.25,0.75,0.13,0.38,0.62,0.88,1.0};

	for (iter=1;iter<=MAXIT;iter++) {
		*its=iter;
		b=a[m];
		err=NR_Cabs(b);
		d=f=NR_Complex(0.0,0.0);
		abx=NR_Cabs(*x);
		for (j=m-1;j>=0;j--) {
			f=NR_Cadd(NR_Cmul(*x,f),d);
			d=NR_Cadd(NR_Cmul(*x,d),b);
			b=NR_Cadd(NR_Cmul(*x,b),a[j]);
			err=NR_Cabs(b)+abx*err;
		}
		err *= EPSS;
		if (NR_Cabs(b) <= err) return;
		g=NR_Cdiv(d,b);
		g2=NR_Cmul(g,g);
		h=NR_Csub(g2,NR_RCmul(2.0,NR_Cdiv(f,b)));
		sq=NR_Csqrt(NR_RCmul((float) (m-1),NR_Csub(NR_RCmul((float) m,h),g2)));
		gp=NR_Cadd(g,sq);
		gm=NR_Csub(g,sq);
		abp=NR_Cabs(gp);
		abm=NR_Cabs(gm);
		if (abp < abm) gp=gm;
		dx=((NR_FMAX(abp,abm) > 0.0 ? NR_Cdiv(NR_Complex((float) m,0.0),gp)
			: NR_RCmul(exp(log(1+abx)),NR_Complex(cos((float)iter),sin((float)iter)))));
		x1=NR_Csub(*x,dx);
		if (x->r == x1.r && x->i == x1.i) return;
		if (iter % MT) *x=x1;
		else *x=NR_Csub(*x,NR_RCmul(frac[iter/MT],dx));
	}
	NR_error("too many iterations in NR_laguer");
	return;
}
#undef EPSS
#undef MR
#undef MT
#undef MAXIT
#undef NRANSI
