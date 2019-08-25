/**
 * @file od_worker.c
 * @brief code for "Solving for optimal defense level" in Method section
 * @author Yuanfei Pan
 * @email yfpan16@fudan.edu.cn
 * @date 2019.08.25
 * @license MIT
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

#include "model.h"
#include "utils.h"
#include "define.h"

int main(int argc, char* argv[])
{
    // timing
    clock_t t = clock();
    
    // define parameters and inputs

                                // Varname      Corresponding parameter name in model
    float x[X_SIZE];            // x        ->  \gamma
    float xa[A_SIZE][X_SIZE];   // xa       ->  \gamma^a
    float xb[B_SIZE][X_SIZE];   // xb       ->  \gamma^b
    float fitness[X_SIZE];      // fitness  ->  W
    float hrng[H_SIZE];         // h        ->  H_0 (`rng` is short for range, `hrng` -> range for parameter H_0)
    float prng[P_SIZE];         // p        ->  p
    float arng[A_SIZE];         // a        ->  a
    float brng[B_SIZE];         // b        ->  b
    float crng[C_SIZE];         // c        ->  c
    float lrng[L_SIZE];         // l        -> \lambda

    // initialize
    // function `void linspace(...)` is defined in model.c
    linspace(x, 0, 1, X_SIZE);
    linspace(hrng, 0, 2, H_SIZE);
    linspace(prng, 0, 1, P_SIZE);
    linspace(arng, 1, 3, A_SIZE);
    linspace(brng, 0, 3, B_SIZE);
    linspace(crng, 0, 1, C_SIZE);
    linspace(lrng, -1, 1, L_SIZE);

    // precompute power(x,a) and power(x,b) to accelerating computation (build lookup table)
    for (int xind = 0; xind < X_SIZE; ++xind)
    {
        for (int aind = 0; aind < A_SIZE; ++aind)
            xa[aind][xind] = pow(x[xind], arng[aind]);
        for (int bind = 0; bind < A_SIZE; ++bind)
            xb[bind][xind] = pow(x[xind], brng[bind]);
    }

    // open output file
    char* output_filename;
    FILE* fp;
    
    output_filename = argv[1];
    fp = fopen(output_filename, "w");

    // split parameter space into pieces (accelerate by running in parallel)
    int shift[] = {0, 0, 0, 0, 0, 0};
    int length[] = {A_SIZE, B_SIZE, H_SIZE, P_SIZE, C_SIZE, L_SIZE};
    int split_on;

    if      (strncmp(argv[2], "h", 1) == 0)
        split_on = 2;
    else if (strncmp(argv[2], "p", 1) == 0)
        split_on = 3;
    else if (strncmp(argv[2], "a", 1) == 0)
        split_on = 0;
    else if (strncmp(argv[2], "b", 1) == 0)
        split_on = 1;
    else if (strncmp(argv[2], "c", 1) == 0)
        split_on = 4;
    else if (strncmp(argv[2], "l", 1) == 0)
        split_on = 5;
    else
        split_on = -1;

    if (split_on != -1)
    {
        shift[split_on] = atoi(argv[3]);
        length[split_on] = atoi(argv[4]);
    }
    
    // start computation
    float param[4];
    int ind_max;
    int progress, total;
    progress = 0;
    total = 1;

    for (int i = 0; i < 6; i++)
        total *= length[i];

    // write the header of output .csv file
    fprintf(fp, "h,p,a,b,c,l\n");
    
    // for each parameter settings
    for (int aind = shift[0]; aind < shift[0]+length[0]; ++aind)
    {
        for (int bind = shift[1]; bind < shift[1]+length[1]; ++bind)
        {
            for (int hind = shift[2]; hind < shift[2]+length[2]; ++hind)
            {
                param[0] = hrng[hind];
                for (int pind = shift[3]; pind < shift[3]+length[3]; ++pind)
                {
                    param[1] = prng[pind];
                    for (int cind = shift[4]; cind < shift[4]+length[4]; ++cind)
                    {
                        param[2] = crng[cind];
                        for (int lind = shift[5]; lind < shift[5]+length[5]; ++lind)
                        {
                            param[3] = lrng[lind];
                            // function `void model(...)` is defined in model.c
                            model(fitness, x, xa[aind], xb[bind], param, X_SIZE);
                            // function `int argmax(...)` is defined in model.c
                            ind_max = argmax(fitness, X_SIZE);

                            if (brng[bind] <= 1)
                            {
                                if (ind_max > 0 && ind_max < X_SIZE-1)
                                {
                                    fprintf(fp, "%.2f,%.2f,%.2f,%.2f,%.2f,%.2f\n", hrng[hind], prng[pind], arng[aind], brng[bind], crng[cind], lrng[lind]);
                                }
                            }
                            //++progress;
                            //printf("%.2f%% %d / %d\r", (float)progress/(float)total*100, progress, total);
                        }
                    }
                }
            }
        }
    }
    fclose(fp);
    printf("\ndone in %d sec.\n", (int)(clock()-t)/CLOCKS_PER_SEC);
    return 0;
}
