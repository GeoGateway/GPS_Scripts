#include <math.h>
#define NRANSI
#include "nr_util.h"

void NR_ks2d1s(float x1[], float y1[], unsigned long n1,
	void (*NR_quadvl)(float, float, float *, float *, float *, float *),
	float *d1, float *prob)
{
	void NR_pearsn(float x[], float y[], unsigned long n, float *r, float *prob,
		float *z);
	float NR_probks(float alam);
	void NR_quadct(float x, float y, float xx[], float yy[], unsigned long nn,
		float *fa, float *fb, float *fc, float *fd);
	unsigned long j;
	float dum,dumm,fa,fb,fc,fd,ga,gb,gc,gd,r1,rr,sqen;

	*d1=0.0;
	for (j=1;j<=n1;j++) {
		NR_quadct(x1[j],y1[j],x1,y1,n1,&fa,&fb,&fc,&fd);
		(*NR_quadvl)(x1[j],y1[j],&ga,&gb,&gc,&gd);
		*d1=NR_FMAX(*d1,fabs(fa-ga));
		*d1=NR_FMAX(*d1,fabs(fb-gb));
		*d1=NR_FMAX(*d1,fabs(fc-gc));
		*d1=NR_FMAX(*d1,fabs(fd-gd));
	}
	NR_pearsn(x1,y1,n1,&r1,&dum,&dumm);
	sqen=sqrt((double)n1);
	rr=sqrt(1.0-r1*r1);
	*prob=NR_probks(*d1*sqen/(1.0+rr*(0.25-0.75/sqen)));
}
#undef NRANSI
