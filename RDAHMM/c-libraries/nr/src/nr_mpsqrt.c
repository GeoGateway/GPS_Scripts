#include <math.h>
#define NRANSI
#include "nr_util.h"
#define MF 3
#define BI (1.0/256)

void NR_mpsqrt(unsigned char w[], unsigned char u[], unsigned char v[], int n,
	int m)
{
	void mplsh(unsigned char u[], int n);
	void mpmov(unsigned char u[], unsigned char v[], int n);
	void NR_mpmul(unsigned char w[], unsigned char u[], unsigned char v[], int n,
		int m);
	void mpneg(unsigned char u[], int n);
	void mpsdv(unsigned char w[], unsigned char u[], int n, int iv, int *ir);
	int i,ir,j,mm;
	float fu,fv;
	unsigned char *r,*s;

	r=NR_cvector(1,n<<1);
	s=NR_cvector(1,n<<1);
	mm=NR_IMIN(m,MF);
	fv=(float) v[mm];
	for (j=mm-1;j>=1;j--) {
		fv *= BI;
		fv += v[j];
	}
	fu=1.0/sqrt(fv);
	for (j=1;j<=n;j++) {
		i=(int) fu;
		u[j]=(unsigned char) i;
		fu=256.0*(fu-i);
	}
	for (;;) {
		NR_mpmul(r,u,u,n,n);
		mplsh(r,n);
		NR_mpmul(s,r,v,n,m);
		mplsh(s,n);
		mpneg(s,n);
		s[1] -= 253;
		mpsdv(s,s,n,2,&ir);
		for (j=2;j<n;j++) {
			if (s[j]) {
				NR_mpmul(r,s,u,n,n);
				mpmov(u,&r[1],n);
				break;
			}
		}
		if (j<n) continue;
		NR_mpmul(r,u,v,n,m);
		mpmov(w,&r[1],n);
		NR_free_cvector(s,1,n<<1);
		NR_free_cvector(r,1,n<<1);
		return;
	}
}
#undef MF
#undef BI
#undef NRANSI
