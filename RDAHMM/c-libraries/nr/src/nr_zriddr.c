#include <math.h>
#define NRANSI
#include "nr_util.h"
#define MAXIT 60
#define UNUSED (-1.11e30)

float NR_zriddr(float (*func)(float), float x1, float x2, float xacc)
{
	int j;
	float ans,fh,fl,fm,fnew,s,xh,xl,xm,xnew;

	fl=(*func)(x1);
	fh=(*func)(x2);
	if ((fl > 0.0 && fh < 0.0) || (fl < 0.0 && fh > 0.0)) {
		xl=x1;
		xh=x2;
		ans=UNUSED;
		for (j=1;j<=MAXIT;j++) {
			xm=0.5*(xl+xh);
			fm=(*func)(xm);
			s=sqrt(fm*fm-fl*fh);
			if (s == 0.0) return ans;
			xnew=xm+(xm-xl)*((fl >= fh ? 1.0 : -1.0)*fm/s);
			if (fabs(xnew-ans) <= xacc) return ans;
			ans=xnew;
			fnew=(*func)(ans);
			if (fnew == 0.0) return ans;
			if (NR_SIGN(fm,fnew) != fm) {
				xl=xm;
				fl=fm;
				xh=ans;
				fh=fnew;
			} else if (NR_SIGN(fl,fnew) != fl) {
				xh=ans;
				fh=fnew;
			} else if (NR_SIGN(fh,fnew) != fh) {
				xl=ans;
				fl=fnew;
			} else NR_error("never get here.");
			if (fabs(xh-xl) <= xacc) return ans;
		}
		NR_error("NR_zriddr exceed maximum iterations");
	}
	else {
		if (fl == 0.0) return x1;
		if (fh == 0.0) return x2;
		NR_error("root must be bracketed in NR_zriddr.");
	}
	return 0.0;
}
#undef MAXIT
#undef UNUSED
#undef NRANSI
