// Code by Apoorav Singh Deo
// Github: https://github.com/apoorav-singh

#include<stdio.h>
#include<stdlib.h>
#include<math.h>
// This module feeds the energy values to the program "phase_shift.c"

// This code can be integrated with the main code if needed



int main()
{
    double E, E_beg=0.05, E_end=0.55, n=50, delta_E;

    delta_E = (E_end - E_beg)/n;

    FILE *fp = fopen("energy.dat","w");

    E = E_beg;
    
    for (int i=0;i<n;i++)
    {
        fprintf(fp, "%lf\n", E);
        E += delta_E;
    }

    fclose(fp);
}
