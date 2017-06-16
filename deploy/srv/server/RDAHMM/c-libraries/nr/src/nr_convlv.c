#define NRANSI
#include "nr_util.h"

void NR_convlv(float data[], unsigned long n, float respns[], unsigned long m,
	int isign, float ans[])
{
	void NR_realft(float data[], unsigned long n, int isign);
	void NR_twofft(float data1[], float data2[], float fft1[], float fft2[],
		unsigned long n);
	unsigned long i,no2;
	float dum,mag2,*fft;

	fft=NR_vector(1,n<<1);
	for (i=1;i<=(m-1)/2;i++)
		respns[n+1-i]=respns[m+1-i];
	for (i=(m+3)/2;i<=n-(m-1)/2;i++)
		respns[i]=0.0;
	NR_twofft(data,respns,fft,ans,n);
	no2=n>>1;
	for (i=2;i<=n+2;i+=2) {
		if (isign == 1) {
			ans[i-1]=(fft[i-1]*(dum=ans[i-1])-fft[i]*ans[i])/no2;
			ans[i]=(fft[i]*dum+fft[i-1]*ans[i])/no2;
		} else if (isign == -1) {
			if ((mag2=NR_SQR(ans[i-1])+NR_SQR(ans[i])) == 0.0)
				NR_error("Deconvolving at response zero in NR_convlv");
			ans[i-1]=(fft[i-1]*(dum=ans[i-1])+fft[i]*ans[i])/mag2/no2;
			ans[i]=(fft[i]*dum-fft[i-1]*ans[i])/mag2/no2;
		} else NR_error("No meaning for isign in NR_convlv");
	}
	ans[2]=ans[n+1];
	NR_realft(ans,n,-1);
	NR_free_vector(fft,1,n<<1);
}
#undef NRANSI
