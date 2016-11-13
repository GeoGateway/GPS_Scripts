#include <math.h>
#define NRANSI
#include "nr_util.h"
#define ALF 1.0e-4
#define TOLX 1.0e-7

void NR_lnsrch(int n, float xold[], float fold, float g[], float p[], float x[],
	float *f, float stpmax, int *check, float (*func)(float []))
{
	int i;
	float a,alam,alam2,alamin,b,disc,f2,fold2,rhs1,rhs2,sNR_lope,sum,temp,
		test,tmplam;

	*check=0;
	for (sum=0.0,i=1;i<=n;i++) sum += p[i]*p[i];
	sum=sqrt(sum);
	if (sum > stpmax)
		for (i=1;i<=n;i++) p[i] *= stpmax/sum;
	for (sNR_lope=0.0,i=1;i<=n;i++)
		sNR_lope += g[i]*p[i];
	test=0.0;
	for (i=1;i<=n;i++) {
		temp=fabs(p[i])/NR_FMAX(fabs(xold[i]),1.0);
		if (temp > test) test=temp;
	}
	alamin=TOLX/test;
	alam=1.0;
	for (;;) {
		for (i=1;i<=n;i++) x[i]=xold[i]+alam*p[i];
		*f=(*func)(x);
		if (alam < alamin) {
			for (i=1;i<=n;i++) x[i]=xold[i];
			*check=1;
			return;
		} else if (*f <= fold+ALF*alam*sNR_lope) return;
		else {
			if (alam == 1.0)
				tmplam = -sNR_lope/(2.0*(*f-fold-sNR_lope));
			else {
				rhs1 = *f-fold-alam*sNR_lope;
				rhs2=f2-fold2-alam2*sNR_lope;
				a=(rhs1/(alam*alam)-rhs2/(alam2*alam2))/(alam-alam2);
				b=(-alam2*rhs1/(alam*alam)+alam*rhs2/(alam2*alam2))/(alam-alam2);
				if (a == 0.0) tmplam = -sNR_lope/(2.0*b);
				else {
					disc=b*b-3.0*a*sNR_lope;
					if (disc<0.0) NR_error("Roundoff problem in NR_lnsrch.");
					else tmplam=(-b+sqrt(disc))/(3.0*a);
				}
				if (tmplam>0.5*alam)
					tmplam=0.5*alam;
			}
		}
		alam2=alam;
		f2 = *f;
		fold2=fold;
		alam=NR_FMAX(tmplam,0.1*alam);
	}
}
#undef ALF
#undef TOLX
#undef NRANSI
