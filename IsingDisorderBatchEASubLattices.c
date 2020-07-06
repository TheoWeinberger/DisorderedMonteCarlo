/////////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////////////
/*                                           Author: 8230W                                                 */
/////////////////////////////////////////////////////////////////////////////////////////////////////////////
/*                                          Code Description                                               */
/////////////////////////////////////////////////////////////////////////////////////////////////////////////
/* Code to run monte carlo simulation of the disordered AFM lattice. Runs the simulation over a specified  */
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
/* Thermodynamics quantities, e, m, sublattice magnetisation m1, EA parameter Q                            */
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

void Measure(double *** Spin_Array, int L, double * Ave_Energy, double * Ave_Spin, double * M_1, double * M_2, double * M_3, double * M_4, double Coupling, int N)
{

    /*Reset measurements*/
    *Ave_Energy = 0.0;

    /*summed parameters for calculations*/
    double M = 0.0;
    double M_1_Total = 0.0;
    int N_1 = 0;
    double M_2_Total = 0.0;
    int N_2 = 0;
    double M_3_Total = 0.0;
    int N_3 = 0;
    double M_4_Total = 0.0;
    int N_4 = 0;
    double E = 0.0;
    int i,j,k;
    for (i = 0; i < L; i++)
    {
        for (j = 0; j < L; j++)
        {
            for (k = 0; k < L; k++)
            {
                if(i + j - k >= 0 && i + j - k < L && i - j + k >= 0 && i - j + k < L && - i + j + k >= 0 && - i + j + k < L && Spin_Array[i+j-k][i-j+k][-i+j+k] != 0.0)
                {
                    M_1_Total += Spin_Array[i+j-k][i-j+k][-i+j+k];
                    N_1++;
                }
                if(i + j - k + 1 >= 0 && i + j - k + 1 < L && i - j + k >= 0 && i - j + k < L && - i + j + k >= 0 && - i + j + k < L && Spin_Array[i+j-k+1][i-j+k][-i+j+k] != 0.0)
                {
                    M_2_Total += Spin_Array[i+j-k+1][i-j+k][-i+j+k];
                    N_2++;
                }
                if(i + j - k >= 0 && i + j - k < L && i - j + k + 1 >= 0 && i - j + k + 1< L && - i + j + k >= 0 && - i + j + k < L && Spin_Array[i+j-k][i-j+k+1][-i+j+k] != 0.0)
                {
                    M_3_Total += Spin_Array[i+j-k][i-j+k+1][-i+j+k];
                    N_3++;
                }
                if(i + j - k >= 0 && i + j - k < L && i - j + k >= 0 && i - j + k < L && - i + j + k + 1 >= 0 && - i + j + k + 1< L && Spin_Array[i+j-k][i-j+k][-i+j+k+1] != 0.0)
                {
                    M_4_Total += Spin_Array[i+j-k][i-j+k][-i+j+k+1];
                    N_4++;
                }
                E += Hamiltonian(Spin_Array, L, Coupling, i, j, k)/2; 
                M += Spin_Array[i][j][k];
            }   
        }
    }
    *M_1 = M_1_Total/N_1;
    *M_2 = M_2_Total/N_2;
    *M_3 = M_3_Total/N_3;
    *M_4 = M_4_Total/N_4;
    *Ave_Spin = M/N;
    *Ave_Energy = E/N;
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

void Edwards_Anderson(double *** Spin_Array, double T, int N, int L, double Coupling, double * Q)
{
    /*reset Q*/
    *Q = 0;

    /*Temporary variable*/
    double Q_Sum = 0;

    /*Counting Parameters*/
    int i,j,k,n;

    /*Array to keep summed spin values in*/
    double *** Sum_Spin_Array = (double ***)malloc(L*sizeof(double**));
    for (i = 0; i < L; i++)
        {
            Sum_Spin_Array[i] = (double **)malloc(L*sizeof(double *));
            for (j = 0; j < L; j++)
            {
                Sum_Spin_Array[i][j] = (double*)malloc(L*sizeof(double));
            }
        }

    /*Reset Sum_Spin_Array*/
    for (i = 0; i < L; i++)
    {
        for (j = 0; j < L; j++)
        {
            for (k = 0; k < L; k++)
            {
                Sum_Spin_Array[i][j][k] = 0;
            }   
        }
    }
    
    /*Run 1000 MonteCarlo steps and sum spins on each site*/
    for(n = 0; n < 1000; n++)
    {
        for (i = 0; i < L; i++)
        {
            for (j = 0; j < L; j++)
            {
                for (k = 0; k < L; k++)
                {
                    Sum_Spin_Array[i][j][k] = Sum_Spin_Array[i][j][k] + Spin_Array[i][j][k];
                }   
            }
        }

        MonteCarlo(Spin_Array, T, 1, N, L, Coupling);

    }

    for (i = 0; i < L; i++)
    {
        for (j = 0; j < L; j++)
        {
            for (k = 0; k < L; k++)
            {
                Q_Sum += fabs(Sum_Spin_Array[i][j][k]);
            }   
        }
    }
    *Q = Q_Sum/(N*1000); 
}

int main()
{
    clock_t begin = clock();
    /*seed random number*/
    //seed_rng();
    srand((unsigned)time(NULL));

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
    int N; /*number of spins in the system*/
    double T; /*Dimensionless temperature T/J */
    double Coupling = -1.0;

    /*MonteCarlo Variables*/
    int Num_Cycles; /*Number of cycles where one cycle includes N spin flips*/

    /*define temperature parameters*/
    double T_Init = 4; /*initial Temperature*/
    double T_Final = 0.6; /*final Temperature chosen to be T_c/2. for AFM = 0.8 for FM = 4.5*/
    double Temp_Scaling = pow(T_Final/T_Init, 1.0/100.0); /*scaling factor for temperature such that there are 100 steps between T_Init and T_Final*/

    /*define MC time parameters*/
    double Step_Scaling = 6.0; /*scaling factor for MC steps chosen according to Tm = 10000e^(x/T) such that Tm at T_final = 1000000*/

    /*Measurements*/
    double M_1 = 0.0;
    double M_2 = 0.0;
    double M_3 = 0.0;
    double M_4 = 0.0;
    double Ave_Spin = 0.0;
    double Ave_Energy = 0.0;
    double Q = 0.0;

    int m = 0; /*counting parameters*/

    /*Title columns*/
    printf("Temperature,Average Energy,Average Spin,M1,M2,M3,M4,Q,N\n");

    /*Reset system*/
    T = T_Init; /*set initial temp*/
    Num_Cycles = round(2*exp(Step_Scaling/T)); /*Set initial number of cycles note additional factor of 2 to ensure equilibration*/

    /*Reset spins to initial configuration to average over a system*/
    N = CSVInitialise(Spin_Array, L);

    MonteCarlo(Spin_Array, 100, 10, N, L, Coupling); /*Run a few runs at a high temperature to scramble the system each time*/

    for (m = 0; m < 100; m++)
    {
        /*Run Simulation*/
        MonteCarlo(Spin_Array, T, Num_Cycles, N, L, Coupling);

        /*Calculate Edwards-Anderson parameter*/
        Edwards_Anderson(Spin_Array, T, N, L, Coupling, &Q);

        /*Calculate Observables*/
        Measure(Spin_Array, L, &Ave_Energy, &Ave_Spin, &M_1, &M_2, &M_3, &M_4, Coupling, N);

        /*Store measurements in data file*/
        printf("%f,%f,%f,%f,%f,%f,%f,%f,%d\n", T, Ave_Energy, Ave_Spin, M_1, M_2, M_3, M_4, Q, N);

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

    clock_t end = clock();
    double time_spent = (double)(end - begin) / CLOCKS_PER_SEC;
    //printf("%f", time_spent);
}