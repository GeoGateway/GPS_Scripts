void NR_rank(unsigned long n, unsigned long indx[], unsigned long iNR_rank[])
{
	unsigned long j;

	for (j=1;j<=n;j++) iNR_rank[indx[j]]=j;
}
