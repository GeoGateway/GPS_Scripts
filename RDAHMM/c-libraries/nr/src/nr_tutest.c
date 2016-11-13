#include <math.h>
#define NRANSI
#include "nr_util.h"

void NR_tutest(float data1[], unsigned long n1, float data2[], unsigned long n2,
	float *t, float *prob)
{
	void NR_avevar(float data[], unsigned long n, float *ave, float *var);
	float NR_betai(float a, float b, float x);
	float var1,var2,df,ave1,ave2;

	NR_avevar(data1,n1,&ave1,&var1);
	NR_avevar(data2,n2,&ave2,&var2);
	*t=(ave1-ave2)/sqrt(var1/n1+var2/n2);
	df=NR_SQR(var1/n1+var2/n2)/(NR_SQR(var1/n1)/(n1-1)+NR_SQR(var2/n2)/(n2-1));
	*prob=NR_betai(0.5*df,0.5,df/(df+NR_SQR(*t)));
}
#undef NRANSI
