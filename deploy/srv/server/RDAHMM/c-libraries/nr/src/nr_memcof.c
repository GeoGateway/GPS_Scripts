#include <math.h>
#define NRANSI
#include "nr_util.h"

void NR_memcof(float data[], int n, int m, float *xms, float d[])
{
	int k,j,i;
	float p=0.0,*wk1,*wk2,*wkm;

	wk1=NR_vector(1,n);
	wk2=NR_vector(1,n);
	wkm=NR_vector(1,m);
	for (j=1;j<=n;j++) p += NR_SQR(data[j]);
	*xms=p/n;
	wk1[1]=data[1];
	wk2[n-1]=data[n];
	for (j=2;j<=n-1;j++) {
		wk1[j]=data[j];
		wk2[j-1]=data[j];
	}
	for (k=1;k<=m;k++) {
		float num=0.0,denom=0.0;
		for (j=1;j<=(n-k);j++) {
			num += wk1[j]*wk2[j];
			denom += NR_SQR(wk1[j])+NR_SQR(wk2[j]);
		}
		d[k]=2.0*num/denom;
		*xms *= (1.0-NR_SQR(d[k]));
		for (i=1;i<=(k-1);i++)
			d[i]=wkm[i]-d[k]*wkm[k-i];
		if (k == m) {
			NR_free_vector(wkm,1,m);
			NR_free_vector(wk2,1,n);
			NR_free_vector(wk1,1,n);
			return;
		}
		for (i=1;i<=k;i++) wkm[i]=d[i];
		for (j=1;j<=(n-k-1);j++) {
			wk1[j] -= wkm[k]*wk2[j];
			wk2[j]=wk2[j+1]-wkm[k]*wk1[j+1];
		}
	}
	NR_error("never get here in NR_memcof.");
}
#undef NRANSI
