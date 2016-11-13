#include <math.h>

float NR_expdev(long *idum)
{
	float NR_ran1(long *idum);
	float dum;

	do
		dum=NR_ran1(idum);
	while (dum == 0.0);
	return -log(dum);
}
