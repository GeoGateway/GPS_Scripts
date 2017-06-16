#define NRANSI
#include "nr_util.h"
#define MF 4
#define BI (1.0/256)

void NR_mpinv(unsigned char u[], unsigned char v[], int n, int m)
{
	void mpmov(unsigned char u[], unsigned char v[], int n);
	void NR_mpmul(unsigned char w[], unsigned char u[], unsigned char v[], int n,
		int m);
	void mpneg(unsigned char u[], int n);
	unsigned char *rr,*s;
	int i,j,maxmn,mm;
	float fu,fv;

	maxmn=NR_IMAX(n,m);
	rr=NR_cvector(1,1+(maxmn<<1));
	s=NR_cvector(1,maxmn);
	mm=NR_IMIN(MF,m);
	fv=(float) v[mm];
	for (j=mm-1;j>=1;j--) {
		fv *= BI;
		fv += v[j];
	}
	fu=1.0/fv;
	for (j=1;j<=n;j++) {
		i=(int) fu;
		u[j]=(unsigned char) i;
		fu=256.0*(fu-i);
	}
	for (;;) {
		NR_mpmul(rr,u,v,n,m);
		mpmov(s,&rr[1],n);
		mpneg(s,n);
		s[1] -= 254;
		NR_mpmul(rr,s,u,n,n);
		mpmov(u,&rr[1],n);
		for (j=2;j<n;j++)
			if (s[j]) break;
		if (j==n) {
			NR_free_cvector(s,1,maxmn);
			NR_free_cvector(rr,1,1+(maxmn<<1));
			return;
		}
	}
}
#undef MF
#undef BI
#undef NRANSI
