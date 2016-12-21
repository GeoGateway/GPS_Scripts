void NR_hpsel(unsigned long m, unsigned long n, float arr[], float heap[])
{
	void NR_sort(unsigned long n, float arr[]);
	void NR_error(char error_text[]);
	unsigned long i,j,k;
	float swap;

	if (m > n/2 || m < 1) NR_error("probable misuse of NR_hpsel");
	for (i=1;i<=m;i++) heap[i]=arr[i];
	NR_sort(m,heap);
	for (i=m+1;i<=n;i++) {
		if (arr[i] > heap[1]) {
			heap[1]=arr[i];
			for (j=1;;) {
				k=j << 1;
				if (k > m) break;
				if (k != m && heap[k] > heap[k+1]) k++;
				if (heap[j] <= heap[k]) break;
				swap=heap[k];
				heap[k]=heap[j];
				heap[j]=swap;
				j=k;
			}
		}
	}
}
