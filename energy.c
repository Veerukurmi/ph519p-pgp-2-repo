#include<stdio.h>
#include<stdlib.h>
#include<math.h>

int main()
{
    double E, E_beg=0.05, E_end=0.5, n=50, delta_E;

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
