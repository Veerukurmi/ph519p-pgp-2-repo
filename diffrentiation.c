// This function runs differentiation 

#include<math.h>
#include<stdlib.h>
#include<stdio.h>

#include "line_counting.c"

int main()
{
    int n, i, j;

    n = read_file_line("ph_sh.txt");
    double x[n], f[n], d_f[n], h;

    FILE *fp = fopen("ph_sh.txt", "r");
    
    for(i=0;i<n;i++)
    {   
        fscanf(fp, "%lf %lf\n", &x[i], &f[i]);
        // printf("%lf %lf\n", x[i], f[i]);
    }

    fclose(fp);

    h = x[2]-x[1];


    // Diffrentiating using Euler's Method
    FILE *fp1 = fopen("time_delay.txt","w");

    for (i=0;i<(n-2);i++)
    {
        d_f[i] = 2*(f[i+1]-f[i])/h;
        fprintf(fp1,"%lf %lf\n", x[i], d_f[i]);
    }

    fclose(fp1);

}
