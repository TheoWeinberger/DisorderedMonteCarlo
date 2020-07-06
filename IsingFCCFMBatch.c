/////////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////////////
/*                                           Author: 8230W                                                 */
/////////////////////////////////////////////////////////////////////////////////////////////////////////////
/*                                          Code Description                                               */
/////////////////////////////////////////////////////////////////////////////////////////////////////////////
/* Code to run monte carlo simulation of the FCC FM lattice. Runs the simulation over a specified         */
/* temperature range with 100 logarithmically spaced steps over range. The system is equilibrated at each  */
/* temperature and thermodynamic variables are measured. Meant to be run on the CSD3 with many iterations  */
/* in parallel                                                                                             */
/////////////////////////////////////////////////////////////////////////////////////////////////////////////
/*                                               Inputs                                                    */
/////////////////////////////////////////////////////////////////////////////////////////////////////////////
/* Specify system size L to be studied, specify temperature range T for the system to be studied over      */
/////////////////////////////////////////////////////////////////////////////////////////////////////////////
/*                                               Outputs                                                   */
/////////////////////////////////////////////////////////////////////////////////////////////////////////////
/* Thermodynamics quantities, e, m                            */
/////////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////////////

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <fcntl.h>
#include <unistd.h>

/*Seed the random number generator.*/
void seed_rng(void)
{
    int fp = open("/dev/urandom", O_RDONLY);
    if (fp == -1) abort();
    unsigned seed;
    unsigned pos = 0;
    while (pos < sizeof(seed)) {
        int amt = read(fp, (char *) &seed + pos, sizeof(seed) - pos);
        if (amt <= 0) abort();
        pos += amt;
    }
    srand(seed);
    close(fp);
}


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

/* function which initialises the spin array randomly assigning them as either up or down*/
void initialise( double *** Spin_Array, int L)
{
    int i,j,k;
    for (i = 0; i < L; i++)
    {
        for (j = 0; j < L; j++)
        {
            for (k = 0; k < L; k++)
            {
                double rand_num = rand_db();
                if(rand_num <= 0.5)
                {
                    Spin_Array[i][j][k] = -1.0; 
                }
                else
                {
                    Spin_Array[i][j][k] = 1.0;
                }

            }
        }
    }
}

void Measure(double *** Spin_Array, int L, double * Ave_Spin, double * Ave_Energy, double Coupling, int N)
{

    /*Reset measurements*/
    *Ave_Energy = 0.0;
    *Ave_Spin = 0.0;

    /*summed parameters for calculations*/
    double M = 0.0;
    double E = 0.0;
    int i,j,k;
    for (i = 0; i < L; i++)
    {
        for (j = 0; j < L; j++)
        {
            for (k = 0; k < L; k++)
            {
                M += Spin_Array[i][j][k];
                E += Hamiltonian(Spin_Array, L, Coupling, i, j, k)/2;
            }   
        }
    }
    *Ave_Spin = M/N;
    *Ave_Energy = E/N;
}

void MonteCarlo(double *** Spin_Array, double T, int Num_Cycles, int N, int L, double Coupling, double * Ave_Spin, double * Ave_Energy)
{
    /*Monte Carlo parameters*/
    double Energy_Diff; /*Energy difference between states*/
    double Boltz; /*Boltzmann Factor*/
    double rand_num; /*Random Number*/
    int i,j,k,a,b; /*Counting Variables*/

    /*Metropolis Algorithm*/
    for (a = 0; a < Num_Cycles; a++)
    {
        for (b = 0; b < N; b++)
        {
            /*Pick random spin*/
            i = rand() % L;
            j = rand() % L;
            k = rand() % L;

            /*Calculate Energy Difference*/
            Energy_Diff = -2.0*Hamiltonian(Spin_Array, L, Coupling, i, j, k);
            
            /*Calculate Bolztmann Factor*/
            Boltz = exp(-Energy_Diff/T);
            
            /*Decide whether to flip spin*/
            rand_num = rand_db(); 
            if (rand_num < Boltz)
            {
                Spin_Array[i][j][k]*=-1;
            }
        }
    }

    /*Calculate Observables*/
    Measure(Spin_Array, L, Ave_Spin, Ave_Energy, Coupling, N);
}


int main()
{
    /*seed random number*/
    seed_rng();

    /*allocate array in memory for a cubic system of size L^3*/
    int L = 20;
    int i,j,k;

    /*array that will be used to run simulations*/    
    double *** Spin_Array = (double ***)malloc(L*sizeof(double**));
    for (i = 0; i < L; i++)
        {
            Spin_Array[i] = (double **)malloc(L*sizeof(double *));
            for (j = 0; j < L; j++)
            {
                Spin_Array[i][j] = (double*)malloc(L*sizeof(double));
            }
        }
    
    /* define system parameters */
    int N = L*L*L; /*number of spins in the system*/
    double T; /*Dimensionless temperature T/J */
    double Coupling = -1.0;

    /*MonteCarlo Variables*/
    int Num_Cycles; /*Number of cycles where one cycle includes N spin flips*/

    /*define temperature parameters*/
    double T_Init = 13; /*initial Temperature*/
    double T_Final = 8; /*final Temperature chosen to be T_c/2. for AFM = 0.8 for FM = 4.5*/
    double Temp_Scaling = pow(T_Final/T_Init, 1.0/100.0); /*scaling factor for temperature such that there are 100 steps between T_Init and T_Final*/

    /*define MC time parameters*/
    double Step_Scaling = 50.0; /*scaling factor for MC steps chosen according to Tm = 10000e^(x/T) such that Tm at T_final = 1000000*/

    /*Measurements*/
    double Ave_Spin = 0.0;
    double Ave_Energy = 0.0;

    int m = 0; /*counting parameters*/

    /*Title columns*/
    printf("Temperature,Average Energy,Average Spin,N\n");

    /*Reset system*/
    T = T_Init; /*set initial temp*/
    Num_Cycles = round(2*exp(Step_Scaling/T)); /*Set initial number of cycles note additional factor of 2 to ensure equilibration*/

    initialise(Spin_Array, L);

    MonteCarlo(Spin_Array, 100, 10, N, L, Coupling, &Ave_Spin, &Ave_Energy); /*Run a few runs at a high temperature to scramble the system each time*/

    for (m = 0; m < 100; m++)
    {
        /*Run Simulation*/
        MonteCarlo(Spin_Array, T, Num_Cycles, N, L, Coupling, &Ave_Spin, &Ave_Energy);

        /*Store measurements in data file*/
        printf("%f,%f,%f,%d\n", T, Ave_Energy, Ave_Spin, N);

        /*Step Temperature*/
        T*=Temp_Scaling;

        /*Calculate number of steps for next cycle*/
        Num_Cycles = exp(Step_Scaling/T);

    }       

    return 0;
}