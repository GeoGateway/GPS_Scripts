extern long idum;

void NR_ranpt(float pt[], float regn[], int n)
{
	float NR_ran1(long *idum);
	int j;

	for (j=1;j<=n;j++)
		pt[j]=regn[j]+(regn[n+j]-regn[j])*NR_ran1(&idum);
}
