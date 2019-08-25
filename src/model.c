/**
 * @file model.c
 * @author Yuanfei Pan
 * @email yfpan16@fudan.edu.cn
 * @date 2019.08.25
 * @license MIT
 */

#include "model.h"

/**
 * @brief Implementation of Equation (8) in the paper. This function calculates fitness values (W) corresponds to each defense level (\gamma) given a certain parameter setting.
 *
 * @param fitness caculated fitness values are stored here
 * @param x defense level \gamma
 * @param xa precomputed \gamma^a
 * @param xb precomputed \gamma^b
 * @param param parameter, ordered in h, p, l, c (H_0, p, \lambda, c)
 * @return void
 */
void model(float *fitness, const float *x, const float *xa, const float *xb, const float *param, int size)
{
    float h, p, c, l;
    float pg, ps;

    h = param[0];
    p = param[1];
    c = param[2];
    l = param[3];

    for (int i = 0; i < size; ++i)
    {
        pg = p * (1 - x[i]);
        ps = (1 - p) * (1 - xa[i] + c * (x[i] - xa[i]));
        fitness[i] = -xb[i] - h * (pg + ps + l * pg * ps);
    }
}
