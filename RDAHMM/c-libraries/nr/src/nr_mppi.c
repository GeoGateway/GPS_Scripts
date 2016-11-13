#include <stdio.h>
#define NRANSI
#include "nr_util.h"
#define IAOFF 48

void NR_mppi(int n)
{
	void NR_mp2dfr(unsigned char a[], unsigned char s[], int n, int *m);
	void mpadd(unsigned char w[], unsigned char u[], unsigned char v[], int n);
	void NR_mpinv(unsigned char u[], unsigned char v[], int n, int m);
	void mplsh(unsigned char u[], int n);
	void mpmov(unsigned char u[], unsigned char v[], int n);
	void NR_mpmul(unsigned char w[], unsigned char u[], unsigned char v[], int n,
		int m);
	void mpsdv(unsigned char w[], unsigned char u[], int n, int iv, int *ir);
	void NR_mpsqrt(unsigned char w[], unsigned char u[], unsigned char v[], int n,
		int m);
	int ir,j,m;
	unsigned char mm,*x,*y,*sx,*sxi,*t,*s,*pi;

	x=NR_cvector(1,n+1);
	y=NR_cvector(1,n<<1);
	sx=NR_cvector(1,n);
	sxi=NR_cvector(1,n);
	t=NR_cvector(1,n<<1);
	s=NR_cvector(1,3*n);
	pi=NR_cvector(1,n+1);
	t[1]=2;
	for (j=2;j<=n;j++) t[j]=0;
	NR_mpsqrt(x,x,t,n,n);
	mpadd(pi,t,x,n);
	mplsh(pi,n);
	NR_mpsqrt(sx,sxi,x,n,n);
	mpmov(y,sx,n);
	for (;;) {
		mpadd(x,sx,sxi,n);
		mpsdv(x,&x[1],n,2,&ir);
		NR_mpsqrt(sx,sxi,x,n,n);
		NR_mpmul(t,y,sx,n,n);
		mpadd(&t[1],&t[1],sxi,n);
		x[1]++;
		y[1]++;
		NR_mpinv(s,y,n,n);
		NR_mpmul(y,&t[2],s,n,n);
		mplsh(y,n);
		NR_mpmul(t,x,s,n,n);
		mm=t[2]-1;
		for (j=3;j<=n;j++) {
			if (t[j] != mm) break;
		}
		m=t[n+1]-mm;
		if (j <= n || m > 1 || m < -1) {
				NR_mpmul(s,pi,&t[1],n,n);
				mpmov(pi,&s[1],n);
				continue;
		}
		printf("pi=\n");
		s[1]=pi[1]+IAOFF;
		s[2]='.';
		m=mm;
		NR_mp2dfr(&pi[1],&s[2],n-1,&m);
		s[m+3]=0;
		printf(" %64s\n",&s[1]);
		NR_free_cvector(pi,1,n+1);
		NR_free_cvector(s,1,3*n);
		NR_free_cvector(t,1,n<<1);
		NR_free_cvector(sxi,1,n);
		NR_free_cvector(sx,1,n);
		NR_free_cvector(y,1,n<<1);
		NR_free_cvector(x,1,n+1);
		return;
	}
}
#undef IAOFF
#undef NRANSI
