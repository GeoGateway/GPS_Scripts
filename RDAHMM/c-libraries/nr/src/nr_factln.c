float NR_factln(int n)
{
	float NR_gammln(float xx);
	void NR_error(char error_text[]);
	static float a[101];

	if (n < 0) NR_error("Negative factorial in routine NR_factln");
	if (n <= 1) return 0.0;
	if (n <= 100) return a[n] ? a[n] : (a[n]=NR_gammln(n+1.0));
	else return NR_gammln(n+1.0);
}
