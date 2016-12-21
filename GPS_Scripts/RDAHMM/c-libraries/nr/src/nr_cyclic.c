#define NRANSI
#include "nr_util.h"

void NR_cyclic(float a[], float b[], float c[], float alpha, float NR_beta,
	float r[], float x[], unsigned long n)
{
	void NR_tridag(float a[], float b[], float c[], float r[], float u[],
		unsigned long n);
	unsigned long i;
	float fact,gamma,*bb,*u,*z;

	if (n <= 2) NR_error("n too small in NR_cyclic");
	bb=NR_vector(1,n);
	u=NR_vector(1,n);
	z=NR_vector(1,n);
	gamma = -b[1];
	bb[1]=b[1]-gamma;
	bb[n]=b[n]-alpha*NR_beta/gamma;
	for (i=2;i<n;i++) bb[i]=b[i];
	NR_tridag(a,bb,c,r,x,n);
	u[1]=gamma;
	u[n]=alpha;
	for (i=2;i<n;i++) u[i]=0.0;
	NR_tridag(a,bb,c,u,z,n);
	fact=(x[1]+NR_beta*x[n]/gamma)/
		(1.0+z[1]+NR_beta*z[n]/gamma);
	for (i=1;i<=n;i++) x[i] -= fact*z[i];
	NR_free_vector(z,1,n);
	NR_free_vector(u,1,n);
	NR_free_vector(bb,1,n);
}
#undef NRANSI
