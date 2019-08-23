#include <math.h>

// an implementation for function `numpy.linspace(...)` in Python package numpy
void linspace(float *arr, float start, float end, int size)
{
    float interval = (end - start) / (size - 1);
    arr[0] = start;
    for (int i = 1; i < size; ++i)
        arr[i] = arr[i - 1] + interval;
}

// an implementation for function `numpy.argmax(...)` in Python package numpy
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

// parameter k is set fixed at k = 1, so there's no k in the following code.
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
