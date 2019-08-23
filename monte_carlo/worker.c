#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

#include "utils.c"

#define A_SIZE 20
#define B_SIZE 20
#define C_SIZE 10
#define L_SIZE 20

int main(int argc, char *argv[])
{
    // timing
    clock_t t = clock();

    // define parameters and inputs

    // Varname      Corresponding parameter name in model
    float x[X_SIZE];          // x        ->  \gamma
    float xa[A_SIZE][X_SIZE]; // xa       ->  \gamma^a
    float xb[B_SIZE][X_SIZE]; // xb       ->  \gamma^b
    float fitness[X_SIZE];    // fitness  ->  W
    float hrng[H_SIZE];       // h        ->  H_0 (`rng` is short for range, `hrng` -> range for parameter H_0)
    float prng[P_SIZE];       // p        ->  p
    float arng[A_SIZE];       // a        ->  a
    float brng[B_SIZE];       // b        ->  b
    float crng[C_SIZE];       // c        ->  c
    float lrng[L_SIZE];       // l        -> \lambda

    float mat[H_SIZE][P_SIZE];
    float diff[H_SIZE][P_SIZE];
    
    float slope[SLOPE_SIZE];
    
    float hsam[SAMPLE_SIZE];
    float psam[SAMPLE_SIZE];
    int dsam[SAMPLE_SIZE];
    int scen[SAMPLE_SIZE];
    
    float result[6];

    // initialize
    // function `void linspace(...)` is defined in model.c
    linspace(x, 0, 1, X_SIZE);
    linspace(hrng, 0, 2, H_SIZE);
    linspace(prng, 0, 1, P_SIZE);
    linspace(arng, 1, 3, A_SIZE);
    linspace(brng, 1, 3, B_SIZE);
    linspace(crng, 0, 1, C_SIZE);
    linspace(lrng, -1, 1, L_SIZE);
    
    linspace(slope, -10, 10, SLOPE_SIZE);

    // precompute power(x,a) and power(x,b) to accelerating computation (build lookup table)
    for (int xind = 0; xind < X_SIZE; ++xind)
    {
        for (int aind = 0; aind < A_SIZE; ++aind)
            xa[aind][xind] = pow(x[xind], arng[aind]);
        for (int bind = 0; bind < A_SIZE; ++bind)
            xb[bind][xind] = pow(x[xind], brng[bind]);
    }

    // open output file
    char *output_filename;
    FILE *fp;

    output_filename = argv[1];
    fp = fopen(output_filename, "w");

    // split parameter space into pieces (accelerate by running in parallel)
    int shift[] = {0, 0, 0, 0};
    int length[] = {A_SIZE, B_SIZE, C_SIZE, L_SIZE};
    int split_on;

    if (strncmp(argv[2], "a", 1) == 0)
        split_on = 0;
    else if (strncmp(argv[2], "b", 1) == 0)
        split_on = 1;
    else if (strncmp(argv[2], "c", 1) == 0)
        split_on = 2;
    else if (strncmp(argv[2], "l", 1) == 0)
        split_on = 3;
    else
        split_on = -1;

    if (split_on != -1)
    {
        shift[split_on] = atoi(argv[3]);
        length[split_on] = atoi(argv[4]);
    }

    // start computation
    int progress, total;
    progress = 0;
    total = 1;

    for (int i = 0; i < 4; i++)
        total *= length[i];

    // write the header of output .csv file
    fprintf(fp, "a,b,c,l,accu_a,slope_a,accu_b,slope_b,accu_c,slope_c\n");

    for (int aind = shift[0]; aind < shift[0] + length[0]; ++aind)
    {
        for (int bind = shift[1]; bind < shift[1] + length[1]; ++bind)
        {
            for (int cind = shift[2]; cind < shift[2] + length[2]; ++cind)
            {
                for (int lind = shift[3]; lind < shift[3] + length[3]; ++lind)
                {
                    optimal_surface(mat, hrng, prng, x, xa[aind], xb[bind], crng[cind], lrng[lind]);
                    differentiate(diff, mat);
                    monte_carlo(hsam, psam, dsam, scen, mat, diff, hrng, prng);
                    eval(result, hsam, psam, dsam, scen, slope);
                    fprintf(fp, "%.2f,%.2f,%.2f,%.2f,%.2f,%.2f,%.2f,%.2f,%.2f,%.2f\n", arng[aind],  brng[bind], crng[cind], lrng[lind], result[0], result[1], result[2], result[3], result[4], result[5]);
                    //++progress;
                    //printf("%.2f%% %d / %d\r", (float)progress/(float)total*100, progress, total);
                }
            }
        }
    }
    fclose(fp);
    printf("\ndone in %d sec.\n", (int)(clock() - t) / CLOCKS_PER_SEC);
    return 0;
}
