/////////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////////////
/*                                           Author: 8230W                                                 */
/////////////////////////////////////////////////////////////////////////////////////////////////////////////
/*                                          Code Description                                               */
/////////////////////////////////////////////////////////////////////////////////////////////////////////////
/* Code produces a lattice and computers the autocorrelation time, defined as the time for the             */
/* autocorrelation function to reach a value of 1/10. Averages over 25 runs after an initial 5 runs for    */
/* system scrambling. Additionally measures number of times each spin has been flipped, spin type and      */
/* and number of neighbours for each spin                                                                  */
/////////////////////////////////////////////////////////////////////////////////////////////////////////////
/*                                               Inputs                                                    */
/////////////////////////////////////////////////////////////////////////////////////////////////////////////
/* Specify system size L to be studied, specify temperature T for the system to be studied at              */
/////////////////////////////////////////////////////////////////////////////////////////////////////////////
/*                                               Outputs                                                   */
/////////////////////////////////////////////////////////////////////////////////////////////////////////////
/* csv file of autocorrelation function against monte carlo step number, csv files containing an array over*/
/* the lattice of normalised flippability, spin type and number of neighbours                              */
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

void MonteCarlo_Spin_Flip(double *** Spin_Array, double *** Spin_Flip_Count, double *** Spin_Attempt_Flip_Count, double T, int Num_Cycles, int N, int L, double Coupling)
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

                Spin_Attempt_Flip_Count[i][j][k] +=1;
            
                /*Decide whether to flip spin*/
                rand_num = rand_db(); 
                if (rand_num < Boltz)
                {
                    Spin_Array[i][j][k]*=-1;
                    Spin_Flip_Count[i][j][k] += 1;
                }

            //step counter variable if the selected lattice point is a spin rather than non-magnetic 
            b++;      
            }
        }
    }
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

    double *** Spin_Flip_Count = (double ***)malloc(L*sizeof(double**));
    for (i = 0; i < L; i++)
        {
            Spin_Flip_Count[i] = (double **)malloc(L*sizeof(double *));
            for (j = 0; j < L; j++)
            {
                Spin_Flip_Count[i][j] = (double*)malloc(L*sizeof(double));
            }
        }

    double *** Spin_Attempt_Flip_Count = (double ***)malloc(L*sizeof(double**));
    for (i = 0; i < L; i++)
        {
            Spin_Attempt_Flip_Count[i] = (double **)malloc(L*sizeof(double *));
            for (j = 0; j < L; j++)
            {
                Spin_Attempt_Flip_Count[i][j] = (double*)malloc(L*sizeof(double));
            }
        }

    /* define system parameters */
    int N ; /*number of spins in the system*/
    double T = 0.5; /*Dimensionless temperature T/J */
    double Coupling = -1.0;
    int Num_Histories = 50;

    /*Measurements*/
    double C_Function; /*correlation function*/
    int t_mc; /*montecarlo steps*/

    /*Assign arrays in which to store measurements over multiple histories*/
    int n; /*counting parameters*/
    int m; /*counting parameters*/

    double ** correlation_array = (double **)malloc(Num_Histories*sizeof(double*));
        for (n = 0; n < Num_Histories; n++)
        {
            correlation_array[n] = (double *)malloc(1500*sizeof(double));
        }

    /*open data files to store data*/
    FILE *cp;
    cp = fopen("FINAL_CorrelationT05.csv","w+");
    FILE *lp;
    lp = fopen("FINAL_LatticeT05.csv","w+");
    FILE *sp;
    sp = fopen("FINAL_LatticeSpinsT05.csv","w+");
    FILE *np;
    np = fopen("FINAL_NearestNeighbourT05.csv", "w+");
    //FILE *tp;
    //tp = fopen("TotalNeighbourT05SpinHist.csv", "w+");

    /*Zero Spin Flip Arrays*/
    for (i = 0; i < L; i++)
    {
        for (j = 0; j < L; j++)
        {
            for (k = 0; k < L; k++)
            {
                Spin_Flip_Count[i][j][k] = 0;
                Spin_Attempt_Flip_Count[i][j][k] = 0;
            }
        }
    }

    /*Title columns*/
    fprintf(cp, "t_mc,C\n");

    /*Title columns*/
    fprintf(lp, "i,j,k,Spin Flips\n");

    /*Title columns*/
    fprintf(sp, "i,j,k,Spin\n");

    /*Title columns*/
    fprintf(np, "i,j,k,Neighbours\n");

    /*Title columns*/
    //fprintf(tp, "i,j,k,Neighbour Spin\n");

    N = CSVInitialise(Spin_Array, L);

    for (n = 0; n < Num_Histories; n++)
    {
        /*Scramble spins at high temp*/
        MonteCarlo(Spin_Array, 100, 10, N, L, Coupling);

        /*Scramble at temp to be analysed*/
        MonteCarlo(Spin_Array, T, 10, N, L, Coupling);
        /*Set Initial Spin Array*/
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

        m = 0; /*reset indexing parameter*/

        C_Function = Correlation(Spin_Array, Spin_Array_Initial, N, L);
        t_mc = 0;

        while (C_Function > 0.1)
        {
            correlation_array[n][m] = C_Function;
            MonteCarlo_Spin_Flip(Spin_Array, Spin_Flip_Count, Spin_Attempt_Flip_Count, T, 1, N, L, Coupling);
            C_Function = Correlation(Spin_Array, Spin_Array_Initial, N, L);
            t_mc = t_mc + 1;
            m++;
            if(t_mc>1500)
            {
                break;
            }
        }
        correlation_array[n][m] = C_Function;
    }

    for (m = 0; m < 1500; m++)
    {
        /*Define measurements to be calculated for this system*/
        double correlation_sum = 0.0;
        double correlation_ave = 0.0;

        /*Sum over histories at a given temperature*/
        for (n = 0; n < Num_Histories; n++)
        {
            correlation_sum += correlation_array[n][m]; 
        }
        
        /*calculate averages at a given temperature*/
        correlation_ave = correlation_sum/Num_Histories;
        fprintf(cp, "%d,%f\n", m, correlation_ave);
    }
    
    /*print to output data file*/
    for (i = 0; i < L; i++)
    {
        for (j = 0; j < L; j++)
        {
            for (k = 0; k < L; k++)
            {
                if (Spin_Attempt_Flip_Count[i][j][k] != 0)
                {
                    Spin_Flip_Count[i][j][k] = Spin_Flip_Count[i][j][k]/Spin_Attempt_Flip_Count[i][j][k];
                }
                else
                {
                    Spin_Flip_Count[i][j][k] = 0.0;                                    
                }
                
                double Num_Neighbours = fabs(Spin_Array [i?(i-1):(L-1)][j][k]) + fabs(Spin_Array[(i+1)%L][j][k]) +
                fabs(Spin_Array[i][j?(j-1):(L-1)][k]) + fabs(Spin_Array[i][(j+1)%L][k]) + fabs(Spin_Array [i][j][k?(k-1):(L-1)]) + fabs(Spin_Array[i][j][(k+1)%L]) +
                fabs(Spin_Array [i?(i-1):(L-1)][(j+1)%L][k]) + fabs(Spin_Array[(i+1)%L][j?(j-1):(L-1)][k]) + fabs(Spin_Array[i?(i-1):(L-1)][j][(k+1)%L]) + fabs(Spin_Array[(i+1)%L][j][k?(k-1):(L-1)]) +
                fabs(Spin_Array[i][j?(j-1):(L-1)][(k+1)%L]) + fabs(Spin_Array[i][(j+1)%L][k?(k-1):(L-1)]);

                //double Neighbour_Spin = fabs(Spin_Array [i?(i-1):(L-1)][j][k] + Spin_Array[(i+1)%L][j][k] +
                //Spin_Array[i][j?(j-1):(L-1)][k] + Spin_Array[i][(j+1)%L][k] + Spin_Array [i][j][k?(k-1):(L-1)] + Spin_Array[i][j][(k+1)%L] +
                //Spin_Array [i?(i-1):(L-1)][(j+1)%L][k] + Spin_Array[(i+1)%L][j?(j-1):(L-1)][k] + Spin_Array[i?(i-1):(L-1)][j][(k+1)%L] + Spin_Array[(i+1)%L][j][k?(k-1):(L-1)] +
                //Spin_Array[i][j?(j-1):(L-1)][(k+1)%L] + Spin_Array[i][(j+1)%L][k?(k-1):(L-1)]);

                fprintf(lp, "%d,%d,%d,%f\n", i, j, k, Spin_Flip_Count[i][j][k]);
                fprintf(sp, "%d,%d,%d,%f\n", i, j, k, Spin_Array[i][j][k]);
                fprintf(np, "%d,%d,%d,%f\n", i, j, k, Num_Neighbours);
                //fprintf(tp, "%d,%d,%d,%f\n", i, j, k, Neighbour_Spin);

            }      
        }
    }
    fclose(cp);
    fclose(lp);
    fclose(sp);
    fclose(np);
    //fclose(tp);
    return 0;
}
