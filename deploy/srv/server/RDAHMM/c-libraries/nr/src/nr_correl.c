#define NRANSI
#include "nr_util.h"

void NR_correl(float data1[], float data2[], unsigned long n, float ans[])
{
	void NR_realft(float data[], unsigned long n, int isign);
	void NR_twofft(float data1[], float data2[], float fft1[], float fft2[],
		unsigned long n);
	unsigned long no2,i;
	float dum,*fft;

	fft=NR_vector(1,n<<1);
	NR_twofft(data1,data2,fft,ans,n);
	no2=n>>1;
	for (i=2;i<=n+2;i+=2) {
		ans[i-1]=(fft[i-1]*(dum=ans[i-1])+fft[i]*ans[i])/no2;
		ans[i]=(fft[i]*dum-fft[i-1]*ans[i])/no2;
	}
	ans[2]=ans[n+1];
	NR_realft(ans,n,-1);
	NR_free_vector(fft,1,n<<1);
}
#undef NRANSI
