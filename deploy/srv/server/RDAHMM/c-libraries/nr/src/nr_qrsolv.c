void NR_qrsolv(float **a, int n, float c[], float d[], float b[])
{
	void NR_rsolv(float **a, int n, float d[], float b[]);
	int i,j;
	float sum,tau;

	for (j=1;j<n;j++) {
		for (sum=0.0,i=j;i<=n;i++) sum += a[i][j]*b[i];
		tau=sum/c[j];
		for (i=j;i<=n;i++) b[i] -= tau*a[i][j];
	}
	NR_rsolv(a,n,d,b);
}
