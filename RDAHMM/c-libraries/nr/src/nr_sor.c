#include <math.h>
#define MAXITS 1000
#define EPS 1.0e-5

void NR_sor(double **a, double **b, double **c, double **d, double **e,
	double **f, double **u, int jmax, double NR_rjac)
{
	void NR_error(char error_text[]);
	int ipass,j,jsw,l,lsw,n;
	double anorm,anormf=0.0,omega=1.0,NR_resid;

	for (j=2;j<jmax;j++)
		for (l=2;l<jmax;l++)
			anormf += fabs(f[j][l]);
	for (n=1;n<=MAXITS;n++) {
		anorm=0.0;
		jsw=1;
		for (ipass=1;ipass<=2;ipass++) {
			lsw=jsw;
			for (j=2;j<jmax;j++) {
				for (l=lsw+1;l<jmax;l+=2) {
					NR_resid=a[j][l]*u[j+1][l]
						+b[j][l]*u[j-1][l]
						+c[j][l]*u[j][l+1]
						+d[j][l]*u[j][l-1]
						+e[j][l]*u[j][l]
						-f[j][l];
					anorm += fabs(NR_resid);
					u[j][l] -= omega*NR_resid/e[j][l];
				}
				lsw=3-lsw;
			}
			jsw=3-jsw;
			omega=(n == 1 && ipass == 1 ? 1.0/(1.0-0.5*NR_rjac*NR_rjac) :
				1.0/(1.0-0.25*NR_rjac*NR_rjac*omega));
		}
		if (anorm < EPS*anormf) return;
	}
	NR_error("MAXITS exceeded");
}
#undef MAXITS
#undef EPS
