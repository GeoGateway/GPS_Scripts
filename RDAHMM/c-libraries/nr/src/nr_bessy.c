float NR_bessy(int n, float x)
{
	float NR_bessy0(float x);
	float NR_bessy1(float x);
	void NR_error(char error_text[]);
	int j;
	float by,bym,byp,tox;

	if (n < 2) NR_error("Index n less than 2 in NR_bessy");
	tox=2.0/x;
	by=NR_bessy1(x);
	bym=NR_bessy0(x);
	for (j=1;j<n;j++) {
		byp=j*tox*by-bym;
		bym=by;
		by=byp;
	}
	return by;
}
