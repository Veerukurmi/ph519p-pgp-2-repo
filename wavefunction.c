// This program prints the wavefunction for the given potential.

// Numerov's Method is used to solve the Schrodinger's Equation.

#include<stdio.h>
#include<stdlib.h>
#include<math.h>

void main()
{
    // h_bar = Planck's constant by 2*pi
    // m_e = Mass of the electron
    double hbar = 1.0, m_e = 1.0, delta_h;
    double x_beg, x_end, E;

    int  grid, i, j, l;
    // l value for the effective potential

    printf("Enter the l value of the effective Potential (Integer only, Usually 0 to 6): ");
    scanf("%d",&l);

    grid = 3000;
    E = 0.5;    // Energy of the system in Atomic Units      
    x_beg = 1e-9;  // Begining of spataial coordinates
    x_end = 30; // End of spatial coordiantes
    delta_h = (x_end-x_beg)/grid;

    double x[grid], u_wave[grid], f[grid], y[grid], V[grid];
    
    // x = position 
    // u_wave = wavefunction 
    // f = Numerov's function 
    // y = Numerov's function
    // V = Potential

    x[0] = x_beg;

    // Loop 1

    for (i=0;i<grid;i++)
    {
        if(5.8<=x[i] && x[i]<=7.7)
        {
            f[i] = (2*m_e/pow(hbar,2))*(-E+(-0.302)+((hbar/(2*m_e))*(l*(l+1)/pow(x[i],2))));
            y[i] = 1 - (pow(delta_h,2))/(12)*f[i];
            x[i+1] = x[i] + delta_h;
            // V[i] = -0.302;
            // printf("%lf  %lf", x(i), V(i));
        }
        else
        {
            f[i] = (2*m_e/pow(hbar,2))*(-E+((hbar/(2*m_e))*(l*(l+1)/pow(x[i],2))));
            y[i] = 1 - (pow(delta_h,2))/(12)*f[i];
            x[i+1] = x[i] + delta_h;
            // V[i] = 0;
            // printf("%lf  %lf", x(i), V(i));

        } 
    }

    u_wave[0] = pow(x[0],l+1);
    u_wave[1] = pow(x[1],l+1);

    for(i=1;i<grid;i++)
    {
        u_wave[i+1] = (u_wave[i]*(12.0-(10.0*y[i]))-y[i-1]*u_wave[i-1])/y[i+1];
    }

    FILE *fp1 = fopen("wave_data.txt", "w");

    for(i=0;i<grid;i++)
    {
        fprintf(fp1, "%lf %lf\n", x[i], u_wave[i]);
    }

    fclose(fp1);
}
