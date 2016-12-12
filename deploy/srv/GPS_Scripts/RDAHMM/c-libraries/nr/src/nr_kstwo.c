#include <math.h>

void NR_kstwo(float data1[], unsigned long n1, float data2[], unsigned long n2,
	float *d, float *prob)
{
	float NR_probks(float alam);
	void NR_sort(unsigned long n, float arr[]);
	unsigned long j1=1,j2=1;
	float d1,d2,dt,en1,en2,en,fn1=0.0,fn2=0.0;

	NR_sort(n1,data1);
	NR_sort(n2,data2);
	en1=n1;
	en2=n2;
	*d=0.0;
	while (j1 <= n1 && j2 <= n2) {
		if ((d1=data1[j1]) <= (d2=data2[j2])) fn1=j1++/en1;
		if (d2 <= d1) fn2=j2++/en2;
		if ((dt=fabs(fn2-fn1)) > *d) *d=dt;
	}
	en=sqrt(en1*en2/(en1+en2));
	*prob=NR_probks((en+0.12+0.11/en)*(*d));
}
