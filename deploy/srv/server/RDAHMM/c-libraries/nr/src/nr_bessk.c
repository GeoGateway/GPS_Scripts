float NR_bessk(int n, float x)
{
	float NR_bessk0(float x);
	float NR_bessk1(float x);
	void NR_error(char error_text[]);
	int j;
	float bk,bkm,bkp,tox;

	if (n < 2) NR_error("Index n less than 2 in NR_bessk");
	tox=2.0/x;
	bkm=NR_bessk0(x);
	bk=NR_bessk1(x);
	for (j=1;j<n;j++) {
		bkp=bkm+j*tox*bk;
		bkm=bk;
		bk=bkp;
	}
	return bk;
}
