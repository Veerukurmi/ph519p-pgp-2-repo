// This piece of code is works as testing ground for new function

// Not a part of main code

#include<stdio.h>
#include<stdlib.h>
#include<math.h>

#include "bessel_1_2.c"

int main()
{
    double bessel1, bessel2, x=1;

    FILE *fp = fopen("test2.txt","w");


    for(int i=0;i<100;i++)
    {
        bessel1 = jl(2,x);
        bessel2 = nl(2,x);
        fprintf(fp, "%lf    %lf\n", x, bessel2);
        x += 0.1;
    }

    fclose(fp);
}