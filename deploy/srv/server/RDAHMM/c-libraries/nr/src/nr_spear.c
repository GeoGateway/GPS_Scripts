#include <math.h>
#define NRANSI
#include "nr_util.h"

void NR_spear(float data1[], float data2[], unsigned long n, float *d, float *zd,
	float *probd, float *rs, float *probrs)
{
	float NR_betai(float a, float b, float x);
	void NR_crank(unsigned long n, float w[], float *s);
	float NR_erfcc(float x);
	void NR_sort2(unsigned long n, float arr[], float brr[]);
	unsigned long j;
	float vard,t,sg,sf,fac,en3n,en,df,aved,*wksp1,*wksp2;

	wksp1=NR_vector(1,n);
	wksp2=NR_vector(1,n);
	for (j=1;j<=n;j++) {
		wksp1[j]=data1[j];
		wksp2[j]=data2[j];
	}
	NR_sort2(n,wksp1,wksp2);
	NR_crank(n,wksp1,&sf);
	NR_sort2(n,wksp2,wksp1);
	NR_crank(n,wksp2,&sg);
	*d=0.0;
	for (j=1;j<=n;j++)
		*d += NR_SQR(wksp1[j]-wksp2[j]);
	en=n;
	en3n=en*en*en-en;
	aved=en3n/6.0-(sf+sg)/12.0;
	fac=(1.0-sf/en3n)*(1.0-sg/en3n);
	vard=((en-1.0)*en*en*NR_SQR(en+1.0)/36.0)*fac;
	*zd=(*d-aved)/sqrt(vard);
	*probd=NR_erfcc(fabs(*zd)/1.4142136);
	*rs=(1.0-(6.0/en3n)*(*d+(sf+sg)/12.0))/sqrt(fac);
	fac=(*rs+1.0)*(1.0-(*rs));
	if (fac > 0.0) {
		t=(*rs)*sqrt((en-2.0)/fac);
		df=en-2.0;
		*probrs=NR_betai(0.5*df,0.5,df/(df+t*t));
	} else
		*probrs=0.0;
	NR_free_vector(wksp2,1,n);
	NR_free_vector(wksp1,1,n);
}
#undef NRANSI
