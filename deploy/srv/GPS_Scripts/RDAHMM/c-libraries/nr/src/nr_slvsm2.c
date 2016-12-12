#include <math.h>

void NR_slvsm2(double **u, double **rhs)
{
	void NR_fill0(double **u, int n);
	double disc,fact,h=0.5;

	NR_fill0(u,3);
	fact=2.0/(h*h);
	disc=sqrt(fact*fact+rhs[2][2]);
	u[2][2] = -rhs[2][2]/(fact+disc);
}
