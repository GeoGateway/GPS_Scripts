#ifndef lint
static char rcsid[] = "$Id: test_rand.c,v 1.2 1996/04/15 22:51:08 agray Exp $";
#endif
/* $Log: test_rand.c,v $
 * Revision 1.2  1996/04/15  22:51:08  agray
 * added test for utScaledRand().
 * ag,granat
 *
 * Revision 1.1  1996/04/15  22:04:28  agray
 * Initial revision
 *
 * */
 
#include <stdlib.h>
#include <stdio.h>
#include "util.h"

/****************************************************************************** MAIN
 Driver of program.
******************************************************************************/
main(int argc, char **argv)
{
  double bound, randomdouble;
  int randomint, seed;
  
int i;

  /*Test utSeedRandomByClock*/
  printf("Generating a random number from the time using time()\n");
  printf("Seeding the random number generator\n");
  utSeedRandomByClock();
  printf("Generating the random number from the seed\n");
  randomint = rand();
  printf("The random number generated is %d\n",randomint);

  printf("\n");

  /*Test utSeedRandomByTimeOfDay*/
  printf("Generating a random number from the time using gettimeofday()\n");  
  printf("Seeding the random number generator\n"); 
  utSeedRandomByTimeOfDay(); 
  printf("Generating the random number from the seed\n");
  randomint = rand();
  printf("The random number generated is %d\n",randomint);

  printf("\n");

  /*Test utSeedRandom*/
  printf("Generating a random number from a user-determined seed\n");
  printf("Enter the seed to be used: ");
  scanf("%d",&seed);
  printf("Generating the random number from the seed\n");
  randomint = rand();
  printf("The random number generated is %d\n",randomint);
  
  printf("\n");

  /*Test utDrand*/
  printf("Generating a random number in a user-defined range\n");
  printf("Enter the upper bound of the range (lower bound is zero): ");
  scanf("%lf",&bound);

  printf("Generating the random number\n");
  randomdouble = utDRand(bound);
  printf("The random number generated is %lf\n",randomdouble);

  /* Test utScaledRand */
  printf("Generating 10 random numbers in the range 0 to 1\n");
  for (i=0;i<10;i++)
    printf("%d. %f\n", i, utScaledRand(0.0,1.0));
}
