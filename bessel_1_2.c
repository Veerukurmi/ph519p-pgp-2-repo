#include<stdio.h>
#include<stdlib.h>
#include<math.h>


/* jl - beware! unstable recurrence for large l */

double jl ( int l, double x )
{
  double cos(), sin();              /* Builtin functions */
  double jm1, j0, jp1;
  int lm;
 
  jm1= cos (x) / x;
  j0 = sin (x) / x;
  
  for ( lm = 0; lm < l; lm++ ) {
    jp1= (2*lm+1)/x*j0-jm1;
    jm1= j0;
    j0 = jp1;
  }
  return j0;
}

/* nl - beware! unstable recurrence for large l */

double nl ( int l, double x )
{
  double cos(), sin();              /* Builtin functions */
  double nm1, n0, np1;
  int lm;
 
  nm1= sin (x) / x;
  n0 =-cos (x) / x;
  
  for ( lm = 0; lm < l; lm++ ) {
    np1= (2*lm+1)/x*n0-nm1;
    nm1= n0;
    n0 = np1;
  }
  return n0;
}