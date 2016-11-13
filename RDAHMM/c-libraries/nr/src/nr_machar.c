#include <math.h>
#define CONV(i) ((float)(i))

void NR_machar(int *iNR_beta, int *it, int *irnd, int *ngrd, int *machep, int *negep,
	int *iexp, int *minexp, int *maxexp, float *eps, float *epsneg,
	float *xmin, float *xmax)
{
	int i,itemp,iz,j,k,mx,nxres;
	float a,b,NR_beta,NR_betah,NR_betain,one,t,temp,temp1,tempa,two,y,z,zero;

	one=CONV(1);
	two=one+one;
	zero=one-one;
	a=one;
	do {
		a += a;
		temp=a+one;
		temp1=temp-a;
	} while (temp1-one == zero);
	b=one;
	do {
		b += b;
		temp=a+b;
		itemp=(int)(temp-a);
	} while (itemp == 0);
	*iNR_beta=itemp;
	NR_beta=CONV(*iNR_beta);
	*it=0;
	b=one;
	do {
		++(*it);
		b *= NR_beta;
		temp=b+one;
		temp1=temp-b;
	} while (temp1-one == zero);
	*irnd=0;
	NR_betah=NR_beta/two;
	temp=a+NR_betah;
	if (temp-a != zero) *irnd=1;
	tempa=a+NR_beta;
	temp=tempa+NR_betah;
	if (*irnd == 0 && temp-tempa != zero) *irnd=2;
	*negep=(*it)+3;
	NR_betain=one/NR_beta;
	a=one;
	for (i=1;i<=(*negep);i++) a *= NR_betain;
	b=a;
	for (;;) {
		temp=one-a;
		if (temp-one != zero) break;
		a *= NR_beta;
		--(*negep);
	}
	*negep = -(*negep);
	*epsneg=a;
	*machep = -(*it)-3;
	a=b;
	for (;;) {
		temp=one+a;
		if (temp-one != zero) break;
		a *= NR_beta;
		++(*machep);
	}
	*eps=a;
	*ngrd=0;
	temp=one+(*eps);
	if (*irnd == 0 && temp*one-one != zero) *ngrd=1;
	i=0;
	k=1;
	z=NR_betain;
	t=one+(*eps);
	nxres=0;
	for (;;) {
		y=z;
		z=y*y;
		a=z*one;
		temp=z*t;
		if (a+a == zero || fabs(z) >= y) break;
		temp1=temp*NR_betain;
		if (temp1*NR_beta == z) break;
		++i;
		k += k;
	}
	if (*iNR_beta != 10) {
		*iexp=i+1;
		mx=k+k;
	} else {
		*iexp=2;
		iz=(*iNR_beta);
		while (k >= iz) {
			iz *= *iNR_beta;
			++(*iexp);
		}
		mx=iz+iz-1;
	}
	for (;;) {
		*xmin=y;
		y *= NR_betain;
		a=y*one;
		temp=y*t;
		if (a+a != zero && fabs(y) < *xmin) {
			++k;
			temp1=temp*NR_betain;
			if (temp1*NR_beta == y && temp != y) {
				nxres=3;
				*xmin=y;
				break;
			}
		}
		else break;
	}
	*minexp = -k;
	if (mx <= k+k-3 && *iNR_beta != 10) {
		mx += mx;
		++(*iexp);
	}
	*maxexp=mx+(*minexp);
	*irnd += nxres;
	if (*irnd >= 2) *maxexp -= 2;
	i=(*maxexp)+(*minexp);
	if (*iNR_beta == 2 && !i) --(*maxexp);
	if (i > 20) --(*maxexp);
	if (a != y) *maxexp -= 2;
	*xmax=one-(*epsneg);
	if ((*xmax)*one != *xmax) *xmax=one-NR_beta*(*epsneg);
	*xmax /= (*xmin*NR_beta*NR_beta*NR_beta);
	i=(*maxexp)+(*minexp)+3;
	for (j=1;j<=i;j++) {
		if (*iNR_beta == 2) *xmax += *xmax;
		else *xmax *= NR_beta;
	}
}
#undef CONV
