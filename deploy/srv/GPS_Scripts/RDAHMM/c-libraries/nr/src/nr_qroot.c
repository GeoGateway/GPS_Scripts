#include <math.h>
#define NRANSI
#include "nr_util.h"
#define ITMAX 20
#define TINY 1.0e-6

void NR_qroot(float p[], int n, float *b, float *c, float eps)
{
	void NR_poldiv(float u[], int n, float v[], int nv, float q[], float r[]);
	int iter;
	float sc,sb,s,rc,rb,r,dv,delc,delb;
	float *q,*qq,*rem;
	float d[3];

	q=NR_vector(0,n);
	qq=NR_vector(0,n);
	rem=NR_vector(0,n);
	d[2]=1.0;
	for (iter=1;iter<=ITMAX;iter++) {
		d[1]=(*b);
		d[0]=(*c);
		NR_poldiv(p,n,d,2,q,rem);
		s=rem[0];
		r=rem[1];
		NR_poldiv(q,(n-1),d,2,qq,rem);
		sb = -(*c)*(rc = -rem[1]);
		rb = -(*b)*rc+(sc = -rem[0]);
		dv=1.0/(sb*rc-sc*rb);
		delb=(r*sc-s*rc)*dv;
		delc=(-r*sb+s*rb)*dv;
		*b += (delb=(r*sc-s*rc)*dv);
		*c += (delc=(-r*sb+s*rb)*dv);
		if ((fabs(delb) <= eps*fabs(*b) || fabs(*b) < TINY)
			&& (fabs(delc) <= eps*fabs(*c) || fabs(*c) < TINY)) {
			NR_free_vector(rem,0,n);
			NR_free_vector(qq,0,n);
			NR_free_vector(q,0,n);
			return;
		}
	}
	NR_error("Too many iterations in routine NR_qroot");
}
#undef ITMAX
#undef TINY
#undef NRANSI
