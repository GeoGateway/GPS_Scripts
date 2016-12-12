#include "nr_complex.h"
#define ONE NR_Complex(1.0,0.0)

extern nr_fcomplex aa,bb,cc,z0,dz;

void NR_hypdrv(float s, float yy[], float dyyds[])
{
	nr_fcomplex z,y[3],dyds[3];

	y[1]=NR_Complex(yy[1],yy[2]);
	y[2]=NR_Complex(yy[3],yy[4]);
	z=NR_Cadd(z0,NR_RCmul(s,dz));
	dyds[1]=NR_Cmul(y[2],dz);
	dyds[2]=NR_Cmul(NR_Csub(NR_Cmul(NR_Cmul(aa,bb),y[1]),NR_Cmul(NR_Csub(cc,
		NR_Cmul(NR_Cadd(NR_Cadd(aa,bb),ONE),z)),y[2])),
		NR_Cdiv(dz,NR_Cmul(z,NR_Csub(ONE,z))));
	dyyds[1]=dyds[1].r;
	dyyds[2]=dyds[1].i;
	dyyds[3]=dyds[2].r;
	dyyds[4]=dyds[2].i;
}
#undef ONE
