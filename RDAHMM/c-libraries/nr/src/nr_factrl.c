#include <math.h>

float NR_factrl(int n)
{
	float NR_gammln(float xx);
	void NR_error(char error_text[]);
	static int ntop=4;
	static float a[33]={1.0,1.0,2.0,6.0,24.0};
	int j;

	if (n < 0) NR_error("Negative factorial in routine NR_factrl");
	if (n > 32) return exp(NR_gammln(n+1.0));
	while (ntop<n) {
		j=ntop++;
		a[ntop]=a[j]*ntop;
	}
	return a[n];
}
