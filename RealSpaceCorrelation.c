/////////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////////////
/*                                           Author: 8230W                                                 */
/////////////////////////////////////////////////////////////////////////////////////////////////////////////
/*                                          Code Description                                               */
/////////////////////////////////////////////////////////////////////////////////////////////////////////////
/* Code to calculate the real space correlation between spins on the disordered lattice                    */
/////////////////////////////////////////////////////////////////////////////////////////////////////////////
/*                                               Inputs                                                    */
/////////////////////////////////////////////////////////////////////////////////////////////////////////////
/* Specify system size L to be studied, specify the temperature T for the system to be studied at          */
/* Must also import the initial lattice to be studied. Can change the basis spin for correlations to be    */
/* measured against                                                                                        */
/////////////////////////////////////////////////////////////////////////////////////////////////////////////
/*                                               Outputs                                                   */
/////////////////////////////////////////////////////////////////////////////////////////////////////////////
/* csv file containing the real space correlation between spins within the lattice                         */
/////////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////////////

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <fcntl.h>
#include <unistd.h>

#define M_PI 3.14159265358979323846

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

int CSVInitialise(double *** Spin_Array, int L)
{
    /*set number of rows in data file*/
    int NUMROWS = L*L*L;

    /*define structure of data to be stored*/
    typedef struct
    {
        int i;
        int j;
        int k;
        float spin;
    } Lattice;

    
    FILE *fp;
    /*set number of rows in structure*/
    Lattice index[NUMROWS];

    /*counting variables*/
    int m,l,a,b,c;

    /*Number of spins in disordered system*/
    int Num_Spins = 0;

    /*Open file containing lattice data*/
    fp = fopen("DisorderLattice.csv", "r");

    for(m = 0; m < NUMROWS; m++) 
    {
        fscanf(fp, " %d,%d,%d,%f", &index[m].i, &index[m].j, &index[m].k, &index[m].spin);
    }


    for (a = 0; a < L; a++)
    {
        for (b = 0; b < L; b++)
        {
            for (c = 0; c < L; c++)
            {
                for (l = 0; l < NUMROWS ; l++)
                {
                    if(a == index[l].i && b == index[l].j && c == index[l].k)
                    {
                        Spin_Array[a][b][c] = index[l].spin;

                        if(index[l].spin != 0)
                        {
                            Num_Spins = Num_Spins + 1;
                        }
                    }
                }
            }
        }
    }


    fclose(fp);


    return Num_Spins;

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
                //Energy = Hamiltonian(Spin_Array, L, Coupling, i, j, k);
                //Energy_Flipped = - Hamiltonian(Spin_Array, L, Coupling, i, j, k);
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


int main()
{
    //clock_t begin = clock();
    /*seed random number*/
    //seed_rng();
    srand((unsigned)time(NULL));

    /*System Parameters*/
    int L = 20;
    double T;
    double Coupling = -1.0;
    int m,n,i,j,k;
    double q_1, q_2, q_3, q_z;

    /*MonteCarlo Variables*/
    int Num_Cycles; /*Number of cycles where one cycle includes N spin flips*/
    int Num_Averages = 20; 

    /*define temperature parameters*/
    double T_Init = 4; /*initial Temperature*/
    double T_Final = 0.6; /*final Temperature chosen to be T_c/2. for AFM = 0.8 for FM = 4.5*/
    double Temp_Scaling = pow(T_Final/T_Init, 1.0/100.0); /*scaling factor for temperature such that there are 100 steps between T_Init and T_Final*/
    double Step_Scaling = 6.0; /*scaling factor for MC steps chosen according to Tm = 10000e^(x/T) such that Tm at T_final = 1000000*/

    /*Reset system*/
    T = T_Init; /*set initial temp*/
    Num_Cycles = round(2*exp(Step_Scaling/T)); /*Set initial number of cycles note additional factor of 2 to ensure equilibration*/


    /*allocate array in memory for a cubic system of size L^3*/  
    double *** Spin_Array = (double ***)malloc(L*sizeof(double**));
    for (i = 0; i < L; i++)
    {
        Spin_Array[i] = (double **)malloc(L*sizeof(double *));
        for (j = 0; j < L; j++)
        {
            Spin_Array[i][j] = (double*)malloc(L*sizeof(double));
        }
    }

    double *** Correlation_Array = (double ***)malloc(L*sizeof(double**));
    for (i = 0; i < L; i++)
    {
        Correlation_Array[i] = (double **)malloc(L*sizeof(double *));
        for (j = 0; j < L; j++)
        {
            Correlation_Array[i][j] = (double*)malloc(L*sizeof(double));
        }
    }


    /*zero correlation array*/
    for (i = 0; i < L; i++)
    {
        for (j = 0; j < L; j++)
        {
            for (k = 0; k < L; k++)
            {
                Correlation_Array[i][j][k] = 0;
            }   
        }
    }

    int N = CSVInitialise(Spin_Array, L); 

    for (m = 0; m < 100; m++)
    {
        /*Run Simulation*/
        MonteCarlo(Spin_Array, T, Num_Cycles, N, L, Coupling);

        /*Step Temperature*/
        T*=Temp_Scaling;

        /*Calculate number of steps for next cycle*/
        Num_Cycles = exp(Step_Scaling/T);

        /*include cutoff for reasonable run time*/
        if(Num_Cycles > 100000)
        {
            Num_Cycles = 100000;
        }
    }   

    /*Average over histories*/
    for (n = 0; n < Num_Averages; n++)
    {
        /*Run Simulation*/
        MonteCarlo(Spin_Array, T, Num_Cycles, N, L, Coupling);

        for (i = 0; i < L; i++)
        {
            for (j = 0; j < L; j++)
            {
                for (k = 0; k < L; k++)
                {
                    Correlation_Array[i][j][k] += Spin_Array[1][0][0]*Spin_Array[i][j][k]; //specific the basis spin for the correlations to be measured against here
                }   
            }
        }        
    }   

    /*title columns*/
    printf("i,j,k,C\n");

    for (i = 0; i < L; i++)
    {
        for (j = 0; j < L; j++)
        {
            for (k = 0; k < L; k++)
            {
                printf("%d,%d,%d,%f\n", i, j, k, Correlation_Array[i][j][k]/Num_Averages);
            }   
        }
    }

    //clock_t end = clock();
    //double time_spent = (double)(end - begin) / CLOCKS_PER_SEC;
    //printf("%f", time_spent);

    return 0;
}