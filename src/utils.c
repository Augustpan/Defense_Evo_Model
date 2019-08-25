/**
 * @file utils.c
 * @brief here I defined some utility functions
 * @author Yuanfei Pan
 * @email yfpan16@fudan.edu.cn
 * @date 2019.08.25
 * @license MIT
 */

#include "utils.h"

/**
 * @brief an implementation for function `numpy.linspace(...)` in Python package numpy
 *
 * @param arr input array
 * @param start smallest element in the desired array
 * @param end largest element in the desired array
 * @param size length of the input array
 *
 * @return void
 */
void linspace(float *arr, float start, float end, int size)
{
    float interval = (end - start) / (size - 1);
    arr[0] = start;
    for (int i = 1; i < size; ++i)
        arr[i] = arr[i - 1] + interval;
}

/**
 * @brief an implementation for function `numpy.argmax(...)` in Python package numpy
 *
 * @param arr input array
 * @param size length of the input array
 *
 * @return int
 * @retval ind_max index of the largest element in the input array
 */
int argmax(const float *arr, int size)
{
    float val_max = arr[0];
    int ind_max = 0;
    
    for (int i = 1; i < size; ++i)
    {
        if (arr[i] > val_max)
        {
            val_max = arr[i];
            ind_max = i;
        }
    }
    return ind_max;
}

/**
 * @brief generate random integer in range [0, maxval)
 *
 * @param maxval just as its name suggests
 *
 * @return int
 */
int randint(int maxval)
{
    return rand() % maxval;
}

/**
 * @brief a function to calculate partial differentrial along the second axis of a 2D matrix
 *
 * @param diff[H_SIZE][P_SIZE] output maxtrix
 * @param mat[H_SIZE][P_SIZE] input maxtrix
 *
 * @return void
 */
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

/**
 * @brief calculate optimal surface, see Figure 1
 *
 * @param mat[H_SIZE][P_SIZE] output maxtrix of optimal surface
 * @param hrng range of parameter H_0
 * @param prng range of parameter p
 * @param x defense level \gamma
 * @param xa \gamma^a, see function void model(...) in model.c
 * @param xb \gamma^b, see function void model(...) in model.c
 * @param c parameter c
 * @param l parameter \lambda
 *
 * @return void
 */
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

/**
 * @brief calculate optimal surface, see Figure 1
 *
 * @param hsam storing the \Delta H_0 value of each simulation
 * @param psam storing the \Delta p value of each simulation
 * @param dsam storing the \Delta \gamma^* (optimal defense level) value of each simulation
 * @param scenario storing the scenario info of each simulation. Scenario A, B, C corespond to 0, 1, 2 repectively. See Figure 2
 * @param mat[H_SIZE][P_SIZE] partial differential with respect to p
 * @param hrng range of parameter H_0
 * @param prng range of parameter p
 *
 * @return void
 */
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

/**
 * @brief a function to evaluate whether the \Delta H_0 - \Delta p surface can be divided by a incline line into two sections to make a good prediction on the change in defense level. see Figure 2
 *
 * @param result storing the accuracy of prediction and the slope of the demarcation line of the best predictive power
 * @param hsam storing the \Delta H_0 value of each simulation
 * @param psam storing the \Delta p value of each simulation
 * @param dsam storing the \Delta \gamma^* (optimal defense level) value of each simulation
 * @param scenario storing the scenario info of each simulation. Scenario A, B, C corespond to 0, 1, 2 repectively. See Figure 2
 * @param slope storing the candidates of slope of demarcation line, to accelerate calculation
 *
 * @return void
 */
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
