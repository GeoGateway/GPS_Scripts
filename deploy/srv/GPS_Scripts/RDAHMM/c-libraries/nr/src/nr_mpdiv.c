#define NRANSI
#include "nr_util.h"
#define MACC 3

void NR_mpdiv(unsigned char q[], unsigned char r[], unsigned char u[],
	unsigned char v[], int n, int m)
{
	void NR_mpinv(unsigned char u[], unsigned char v[], int n, int m);
	void mpmov(unsigned char u[], unsigned char v[], int n);
	void NR_mpmul(unsigned char w[], unsigned char u[], unsigned char v[], int n,
		int m);
	void mpsub(int *is, unsigned char w[], unsigned char u[], unsigned char v[],
		int n);
	int is;
	unsigned char *rr,*s;

	rr=NR_cvector(1,(n+MACC)<<1);
	s=NR_cvector(1,n+MACC);
	NR_mpinv(s,v,n-m+MACC,m);
	NR_mpmul(rr,s,u,n-m+MACC,n);
	mpmov(q,&rr[1],n-m+1);
	NR_mpmul(rr,q,v,n-m+1,m);
	mpsub(&is,&rr[1],u,&rr[1],n);
	if (is) NR_error("MACC too small in NR_mpdiv");
	mpmov(r,&rr[n-m+1],m);
	NR_free_cvector(s,1,n+MACC);
	NR_free_cvector(rr,1,(n+MACC)<<1);
}
#undef MACC
#undef NRANSI
