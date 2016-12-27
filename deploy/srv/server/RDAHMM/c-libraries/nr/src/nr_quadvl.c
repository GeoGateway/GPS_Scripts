#define NRANSI
#include "nr_util.h"

void NR_quadvl(float x, float y, float *fa, float *fb, float *fc, float *fd)
{
	float qa,qb,qc,qd;

	qa=NR_FMIN(2.0,NR_FMAX(0.0,1.0-x));
	qb=NR_FMIN(2.0,NR_FMAX(0.0,1.0-y));
	qc=NR_FMIN(2.0,NR_FMAX(0.0,x+1.0));
	qd=NR_FMIN(2.0,NR_FMAX(0.0,y+1.0));
	*fa=0.25*qa*qb;
	*fb=0.25*qb*qc;
	*fc=0.25*qc*qd;
	*fd=0.25*qd*qa;
}
#undef NRANSI
