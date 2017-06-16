void NR_slvsml(double **u, double **rhs)
{
	void NR_fill0(double **u, int n);
	double h=0.5;

	NR_fill0(u,3);
	u[2][2] = -h*h*rhs[2][2]/4.0;
}
