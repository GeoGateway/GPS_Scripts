#include <math.h>
#define NRANSI
#include "nr_util.h"
#define M 64
#define NDFT 1024
#define MPOL 6
#define TWOPI (2.0*3.14159265)

void NR_dftint(float (*func)(float), float a, float b, float w, float *cosint,
	float *sinint)
{
	void NR_dftcor(float w, float delta, float a, float b, float endpts[],
		float *corre, float *corim, float *corfac);
	void NR_polint(float xa[], float ya[], int n, float x, float *y, float *dy);
	void NR_realft(float data[], unsigned long n, int isign);
	static int init=0;
	int j,nn;
	static float aold = -1.e30,bold = -1.e30,delta,(*funcold)(float);
	static float data[NDFT+1],endpts[9];
	float c,cdft,cerr,corfac,corim,corre,en,s;
	float sdft,serr,*cpol,*spol,*xpol;

	cpol=NR_vector(1,MPOL);
	spol=NR_vector(1,MPOL);
	xpol=NR_vector(1,MPOL);
	if (init != 1 || a != aold || b != bold || func != funcold) {
		init=1;
		aold=a;
		bold=b;
		funcold=func;
		delta=(b-a)/M;
		for (j=1;j<=M+1;j++)
			data[j]=(*func)(a+(j-1)*delta);
		for (j=M+2;j<=NDFT;j++)
			data[j]=0.0;
		for (j=1;j<=4;j++) {
			endpts[j]=data[j];
			endpts[j+4]=data[M-3+j];
		}
		NR_realft(data,NDFT,1);
		data[2]=0.0;
	}
	en=w*delta*NDFT/TWOPI+1.0;
	nn=NR_IMIN(NR_IMAX((int)(en-0.5*MPOL+1.0),1),NDFT/2-MPOL+1);
	for (j=1;j<=MPOL;j++,nn++) {
		cpol[j]=data[2*nn-1];
		spol[j]=data[2*nn];
		xpol[j]=nn;
	}
	NR_polint(xpol,cpol,MPOL,en,&cdft,&cerr);
	NR_polint(xpol,spol,MPOL,en,&sdft,&serr);
	NR_dftcor(w,delta,a,b,endpts,&corre,&corim,&corfac);
	cdft *= corfac;
	sdft *= corfac;
	cdft += corre;
	sdft += corim;
	c=delta*cos(w*a);
	s=delta*sin(w*a);
	*cosint=c*cdft-s*sdft;
	*sinint=s*cdft+c*sdft;
	NR_free_vector(cpol,1,MPOL);
	NR_free_vector(spol,1,MPOL);
	NR_free_vector(xpol,1,MPOL);
}
#undef M
#undef NDFT
#undef MPOL
#undef TWOPI
#undef NRANSI
