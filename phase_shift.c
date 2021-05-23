// Phase-Shift 


#include<stdio.h>
#include<stdlib.h>
#include<math.h>

#include "line_counting.c"
#include "bessel_1_2.c"

void main()
{
    // h_bar = Planck's constant by 2*pi
    // m_e = Mass of the electron
    double hbar = 1.0, m_e = 1.0;
    double x_beg, x_end;

    int  grid, i, j, l, n_energy;
    // l value for the effective potential

    printf("Enter the l value of the effective Potential (Integer only, Usually 0 to 6): ");
    scanf("%d",&l);

    grid = 5000;
    x_beg = 1.0;  // Begining of spataial coordinates
    // x_end = 30; // End of spatial coordiantes
    // delta_h = (x_end-x_beg)/grid;



    // n_energy will store total Energy data points

    n_energy = read_file_line("energy.dat"); // 50;

    // double delta_E = (0.55-0.05)/50;


    printf("%d\n",n_energy);

    double delta_h[n_energy];
    
    // E[] is energy array

    double E[n_energy];

    FILE *fp = fopen("energy.dat", "r");

    for (j=0;j<n_energy;j++)
    {
        // E[j] = 0.05 + j*delta_E;
        fscanf(fp, "%lf", &E[j]);
        //printf("%d", j);
    }

    //fclose(fp);

    double ind_inc[n_energy], tan_del[n_energy], ph_sh[n_energy];
    double k[n_energy], K[n_energy], inc[n_energy];

    int r1, r2; 

    // r1 is the initial indices
    // r2 is final indices after addition of the lambda/2

    // x[] is the displacement array
    double x[grid][n_energy];

    double u_wave[grid][n_energy], f[grid][n_energy], y[grid][n_energy];

    //double  V[grid][n_energy];
    
    // x = position 
    // u_wave = wavefunction 
    // f = Numerov's function 
    // y = Numerov's function
    // V = Potential

    // Loop 1

    for (i=0;i<n_energy;i++)
    {
        x[0][i] = x_beg;
    }

    for (j=0;j<n_energy;j++)
    {   
        k[j] = pow((2*m_e*E[j]),0.5);
        inc[j] = (M_PI)/k[j];
        x_end = 20 + inc[j];
        delta_h[j] = (x_end-x_beg)/grid;
        ind_inc[j] = inc[j]/delta_h[j];

        for (i=0;i<grid;i++)
        {
            if(5.8<=x[i][j] && x[i][j]<=7.7)
            {
                f[i][j] = (2*m_e/pow(hbar,2))*(-E[j]+(-0.302)+((hbar/(2*m_e))*(l*(l+1)/pow(x[i][j],2))));
                y[i][j] = 1 - (pow(delta_h[j],2))/(12)*f[i][j];
                x[i+1][j] = x[i][j] + delta_h[j];
                // V[i][j] = -0.302;
                // printf("%lf\n", x[i][j]);
            }
            else
            {
                f[i][j] = (2*m_e/pow(hbar,2))*(-E[j]+((hbar/(2*m_e))*(l*(l+1)/pow(x[i][j],2))));
                y[i][j] = 1 - (pow(delta_h[j],2))/(12)*f[i][j];
                x[i+1][j] = x[i][j] + delta_h[j];
                // V[i][j] = 0;
                // printf("%lf\n", x[i][j]);
            } 
        }
    }

    for (j=0; j<n_energy;j++)
    {
        u_wave[0][j] = pow(x[0][j],l+1);
        u_wave[1][j] = pow(x[1][j],l+1);
    }
    

    for (j=0;j<n_energy;j++)
    {
        for (i=1;i<grid;i++)
        {
            u_wave[i+1][j] = (u_wave[i][j]*(12.0-(10.0*y[i][j]))-y[i-1][j]*u_wave[i-1][j])/y[i+1][j];
        }
    }

    for (i=0;i<n_energy;i++)
     {
         x[0][i] = x_beg;
     }

    FILE *fp1 = fopen("wave_data.txt", "w");

    // Redeclaring the value of x[] as earlier it was giving garbage value

    for(i=0;i<grid;i++)
    {
        fprintf(fp1, "%lf %0.6lf\n",x[i][45], u_wave[i][45]);
    }

    fclose(fp1);

    // Calculation of the phase shift

    // r1 = 10;
    int z=1;
    FILE *fp2 = fopen("ph_sh.txt", "w");

    for (j=0;j<n_energy;j++)
    {
        for (i=0;i<grid;i++)
        {
            // printf("%d\n", z++);
            if(7.0<=x[i][j] && x[i][j]<=7.05)
            {
                r1 = i;
                // printf("%d\n", r1); 
                break;
             }
         }
        r2 = r1 + (int)ind_inc[j];
        // printf("%d\n", r2);
        K[j] = (x[r1][j]*u_wave[r2][j])/(x[r2][j]*u_wave[r1][j]);
        tan_del[j] = (K[j]*jl(l,k[j]*x[r1][j]) - jl(l,k[j]*x[r2][j])/(K[j]*nl(l,k[j]*x[r1][j]) - nl(l,k[j]*x[r2][j])));
        ph_sh[j] = atan(tan_del[j]);
        if (ph_sh[j]<0)
        {
            ph_sh[j] = ph_sh[j] + M_PI;
        }
        fprintf(fp2, "%lf %lf\n", E[j], ph_sh[j]);
    }
    fclose(fp2);
}
