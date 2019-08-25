/**
 * @file utils.h
 * @brief here I defined some utility functions
 * @author Yuanfei Pan
 * @email yfpan16@fudan.edu.cn
 * @date 2019.08.25
 * @license MIT
 */

#ifndef utils_h
#define utils_h

#include <stdio.h>
#include <stdlib.h>
#include <time.h>

#include "model.h"
#include "define.h"

void linspace(float *arr, float start, float end, int size);
int argmax(const float *arr, int size);
int randint(int maxval);
void differentiate(float diff[H_SIZE][P_SIZE], const float mat[H_SIZE][P_SIZE]);
void optimal_surface(float mat[H_SIZE][P_SIZE], const float* hrng, const float* prng, const float* x, const float* xa, const float* xb, float c, float l);
void monte_carlo(float* hsam, float* psam, int* dsam, int* scenario, const float mat[H_SIZE][P_SIZE], const float diff[H_SIZE][P_SIZE], const float* hrng, const float* prng);
void eval(float* result, const float* hsam, const float* psam, const int* dsam, const int* scenario, const float* slope);

#endif /* utils_h */
