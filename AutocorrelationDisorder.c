/////////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////////////
/*                                           Author: 8230W                                                 */
/////////////////////////////////////////////////////////////////////////////////////////////////////////////
/*                                          Code Description                                               */
/////////////////////////////////////////////////////////////////////////////////////////////////////////////
/* Code produces a lattice and computers the autocorrelation time, defined as the time for the             */
/* autocorrelation function to reach a value of 1/10. Averages over 25 runs after an initial 5 runs for    */
/* system scrambling. Modified to produce a disordered lattice allowing for study of the disordered system */
/////////////////////////////////////////////////////////////////////////////////////////////////////////////
/*                                               Inputs                                                    */
/////////////////////////////////////////////////////////////////////////////////////////////////////////////
/* Specify system size L to be studied, specify temperature T for the system to be studied at              */
/////////////////////////////////////////////////////////////////////////////////////////////////////////////
/*                                               Outputs                                                   */
/////////////////////////////////////////////////////////////////////////////////////////////////////////////
/* csv file of autocorrelation function against monte carlo step number                                    */
/////////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////////////

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>


/* function to produce a random float between 0 and 1 to be used by Metropolis Algorithm
and determining initial spin configuration*/
double rand_db()
{
    double rand_num = ((double)rand()/(double)(RAND_MAX)) * 1.0;
    return rand_num;
}

/* function to compute the energy of a spin. Energy is computed as the sum
of interactions of nearest neighbours in a Ising spin model with H = -J SUM(s_i*s_j)*/
double Hamiltonian(double *** Spin_Array, int L, double Coupling, int i, int j, int k)
{
    double Energy;
    Energy = - Coupling*(Spin_Array[i][j][k])*(Spin_Array [i?(i-1):(L-1)][j][k] + Spin_Array[(i+1)%L][j][k] +
    Spin_Array[i][j?(j-1):(L-1)][k] + Spin_Array[i][(j+1)%L][k] + Spin_Array [i][j][k?(k-1):(L-1)] + Spin_Array[i][j][(k+1)%L] +
    Spin_Array [i?(i-1):(L-1)][(j+1)%L][k] + Spin_Array[(i+1)%L][j?(j-1):(L-1)][k] + Spin_Array[i?(i-1):(L-1)][j][(k+1)%L] + Spin_Array[(i+1)%L][j][k?(k-1):(L-1)] +
    Spin_Array[i][j?(j-1):(L-1)][(k+1)%L] + Spin_Array[i][(j+1)%L][k?(k-1):(L-1)]);
    return Energy;
}

/* function to compute the Interaction energy. Energy is computed as the sum
of interactions of nearest neighbours in a Ising spin model with H = -J SUM(s_i*s_j)
note only 6 of the 12 nearest neighbours are included to avoid double counting*/
double Interaction(double *** Spin_Array, int L, double Coupling, int i, int j, int k)
{
    double Energy;
    Energy = - Coupling*(Spin_Array[i][j][k])*(Spin_Array[(i+1)%L][j][k] + Spin_Array[i][(j+1)%L][k] + Spin_Array[i][j][(k+1)%L] +
    Spin_Array [i?(i-1):(L-1)][(j+1)%L][k] + Spin_Array[i?(i-1):(L-1)][j][(k+1)%L] + Spin_Array[i?(i-1):(L-1)][j][k]);
    return Energy;
}

/* function which initialises the spin array randomly assigning them as either up or down*/
void initialise( double *** Spin_Array, int L, int * N)
{
    int i,j,k;
    int Num_Spins = 0;
    *N = 0;
    for (i = 0; i < L; i++)
    {
        for (j = 0; j < L; j++)
        {
            for (k = 0; k < L; k++)
            {
                double rand_num = rand_db();
                if(rand_num <= 0.25)
                {
                    Spin_Array[i][j][k] = -1.0; 
                    Num_Spins = Num_Spins + 1;
                }
                else if(rand_num > 0.25 && rand_num <= 0.5)
                {
                    Spin_Array[i][j][k] = 1.0;
                    Num_Spins = Num_Spins + 1;
                }
                else
                {
                    Spin_Array[i][j][k] = 0.0;
                }

            }
        }
    }
    *N = Num_Spins;
}

void MonteCarlo(double *** Spin_Array, double T, int Num_Cycles, int N, int L, double Coupling)
{
    /*Monte Carlo parameters*/
    double Energy_Diff; /*Energy difference between states*/
    double Boltz; /*Boltzmann Factor*/
    double rand_num; /*Random Number*/
    int i,j,k,a,b; /*Counting Variables*/

    /*Metropolis Algorithm*/
    for (a = 0; a < Num_Cycles; a++)
    {
        b = 0;
        while (b < N)
        {
            /*Pick random spin*/
            i = rand() % L;
            j = rand() % L;
            k = rand() % L;

            if(Spin_Array[i][j][k] != 0.0)
            {
                /*Calculate Energy Difference*/
                Energy_Diff = - 2.0*Hamiltonian(Spin_Array, L, Coupling, i, j, k);
            
                /*Calculate Bolztmann Factor*/
                Boltz = exp(-Energy_Diff/T);
            
                /*Decide whether to flip spin*/
                rand_num = rand_db(); 
                if (rand_num < Boltz)
                {
                    Spin_Array[i][j][k]*=-1;
                }

            //step counter variable if the selected lattice point is a spin rather than non-magnetic 
            b++;      
            }
        }
    }
}

double Correlation(double *** Spin_Array, double *** Spin_Array_Initial, int N, int L)
{
    double correlation = 0;
    int i,j,k;
    for (i = 0; i < L; i++)
    {
        for (j = 0; j < L; j++)
        {
            for (k = 0; k < L; k++)
            {
                correlation += Spin_Array_Initial[i][j][k]*Spin_Array[i][j][k];
            }
        }
    }
    correlation = correlation/N;
    return correlation;
}


int main()
{
    /*seed random number generator*/
    srand((unsigned)time(NULL));

    /*allocate array in memory for a cubic system of size L^3*/
    int L = 20;
    int i,j,k;
    double *** Spin_Array = (double ***)malloc(L*sizeof(double**));
    for (i = 0; i < L; i++)
        {
            Spin_Array[i] = (double **)malloc(L*sizeof(double *));
            for (j = 0; j < L; j++)
            {
                Spin_Array[i][j] = (double*)malloc(L*sizeof(double));
            }
        }
    
    double *** Spin_Array_Initial = (double ***)malloc(L*sizeof(double**));
    for (i = 0; i < L; i++)
        {
            Spin_Array_Initial[i] = (double **)malloc(L*sizeof(double *));
            for (j = 0; j < L; j++)
            {
                Spin_Array_Initial[i][j] = (double*)malloc(L*sizeof(double));
            }
        }

    /* define system parameters */
    int N; /*number of spins in the system*/
    double T; /*Dimensionless temperature T/J */
    double Coupling = -1.0;

    /*define temperature parameters*/
    double T_Init = 3; /*initial Temperature*/
    double T_Final = 0.5; /*final Temperature chosen to be T_c/2. for AFM = 0.8 for FM = 4.5*/

    /*Measurements*/
    double C_Function; /*correlation function*/
    int t_mc; /*montecarlo steps*/

    /*Assign arrays in which to store measurements over multiple histories*/
    int n = 0; /*counting parameters*/
    int Num_Runs = 20; /*Number of runs to average over +5*/

    /*create and open ouput data file*/
    FILE *fp;
    fp = fopen("Correlationdata.csv","w+");

    /*Title columns*/
    fprintf(fp, "T_mc,T\n");

    T = T_Init;
    while (T >= T_Final)
    {
        initialise(Spin_Array, L, &N);
        printf("%d\n", N);
        MonteCarlo(Spin_Array, T, 10, N, L, Coupling);
        double t_mc_ave = 0;
        for (n = 0; n < Num_Runs; n++)
        { 
            for (i = 0; i < L; i++)
            {
                for (j = 0; j < L; j++)
                {
                    for (k = 0; k < L; k++)
                    {
                        Spin_Array_Initial[i][j][k] = Spin_Array[i][j][k];
                    }
                }
            }
        
            C_Function = Correlation(Spin_Array, Spin_Array_Initial, N, L);
            t_mc = 0;

            while (C_Function > 0.1)
            {
                MonteCarlo(Spin_Array, T, 1, N, L, Coupling);
                C_Function = Correlation(Spin_Array, Spin_Array_Initial, N, L);
                //printf("correlation: %f\n", correlation);
                t_mc = t_mc + 1;
                if(t_mc>10000)
                {
                    break;
                }
            }

            if(n >= 5)
            {
                t_mc_ave = t_mc_ave + t_mc;
            }
        } 
        t_mc_ave = t_mc_ave/25;
        fprintf(fp, "%f, %f\n", t_mc_ave, T);
        printf("%f, %f\n", t_mc_ave, T);
        T = T - 0.1;  
    }
    fclose(fp);
    return 0;
}
