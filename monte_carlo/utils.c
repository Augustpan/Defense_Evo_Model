#include <stdio.h>
#include <stdlib.h>
#include <time.h>

#include "model.c"

#define bool int
#define false 0
#define true 1

#define H_SIZE 20
#define P_SIZE 10
#define X_SIZE 1000
#define SAMPLE_SIZE 10000
#define SLOPE_SIZE 200

// generate random int in range [0, m)
int randint(int maxval)
{
    return rand() % maxval;
}

void differentiate(float diff[H_SIZE][P_SIZE], const float mat[H_SIZE][P_SIZE])
{
    for (int hind = 0; hind < H_SIZE; ++hind)
    {
        for (int pind = 1; pind < P_SIZE; ++pind)
        {
            diff[hind][pind] = mat[hind][pind] - mat[hind][pind-1];
        }
        diff[hind][0] = diff[hind][1];
    }
}

void optimal_surface(float mat[H_SIZE][P_SIZE], const float* hrng, const float* prng, const float* x, const float* xa, const float* xb, float c, float l)
{
    float fitness[X_SIZE];
    float param[4];
    
    param[2] = c;
    param[3] = l;
    
    for (int hind = 0; hind < H_SIZE; ++hind)
    {
        param[0] = hrng[hind];
        for (int pind = 0; pind < P_SIZE; ++pind)
        {
            param[1] = prng[pind];
            model(fitness, x, xa, xb, param, X_SIZE);
            mat[hind][pind] = x[argmax(fitness, X_SIZE)];
        }
    }
}

void monte_carlo(float* hsam, float* psam, int* dsam, int* scenario, const float mat[H_SIZE][P_SIZE], const float diff[H_SIZE][P_SIZE], const float* hrng, const float* prng)
{
    float delta_defense;
    int h1, h2, p1, p2;
    
    for (int i = 0; i < SAMPLE_SIZE; ++i)
    {
        h1 = randint(H_SIZE);
        h2 = randint(H_SIZE);
        p1 = randint(P_SIZE);
        p2 = randint(P_SIZE);

        if (diff[h1][p1] > 0 && diff[h2][p2] > 0)
            scenario[i] = 2;
        else if (diff[h1][p1] <= 0 && diff[h2][p2] <= 0)
            scenario[i] = 0;
        else
            scenario[i] = 1;
        
        hsam[i] = hrng[h2] - hrng[h1];
        psam[i] = prng[p2] - prng[p1];
        
        delta_defense = mat[h2][p2] - mat[h1][p1];
        if (delta_defense > 0)
            dsam[i] = 1;     // defense level increased
        else if (delta_defense < 0)
            dsam[i] = -1;    // defense level decreased
        else
            dsam[i] = 0;     // defense level unchanged
    }
}

void eval(float* result, const float* hsam, const float* psam, const int* dsam, const int* scenario, const float* slope)
{
    float predicter;
    bool is_correct;
    
    int total_count[] = {0, 0, 0};
    int cor_count[] = {0, 0, 0};
    int max_cor[] = {0, 0, 0};
    float max_slope[3];
    
    for (int i = 0; i < SLOPE_SIZE; ++i)
    {
        for (int k = 0; k < 3; k++)
        {
            cor_count[k] = 0;
            total_count[k] = 0;
        }

        for (int j = 0; j < SAMPLE_SIZE; ++j)
        {
            if (dsam[j] == 0)
                continue;
            
            predicter = slope[i] * psam[j];
            
            if (predicter > hsam[j] && dsam[j] == -1)
                is_correct = true;
            else if (predicter < hsam[j] && dsam[j] == 1)
                is_correct = true;
            else
                is_correct = false;
            
            if (is_correct)
                ++cor_count[scenario[j]];
            
            ++total_count[scenario[j]];
        }
        
        for (int k = 0; k < 3; k++)
            if (cor_count[k] > max_cor[k])
            {
                max_cor[k] = cor_count[k];
                max_slope[k] = slope[i];
            }
    }
    
    result[0] = (float)max_cor[0] / (float)total_count[0];
    result[1] = max_slope[0];
    result[2] = (float)max_cor[1] / (float)total_count[1];
    result[3] = max_slope[1];
    result[4] = (float)max_cor[2] / (float)total_count[2];
    result[5] = max_slope[2];
}
