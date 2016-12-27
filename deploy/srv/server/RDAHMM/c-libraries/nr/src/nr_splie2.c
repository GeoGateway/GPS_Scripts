void NR_splie2(float x1a[], float x2a[], float **ya, int m, int n, float **y2a)
{
	void NR_spline(float x[], float y[], int n, float yp1, float ypn, float y2[]);
	int j;

	for (j=1;j<=m;j++)
		NR_spline(x2a,ya[j],n,1.0e30,1.0e30,y2a[j]);
}
