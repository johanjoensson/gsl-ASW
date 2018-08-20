#include "utils.h"
#include <iostream>

bool operator==(const lm &a, const lm &b)
{
    return (a.l == b.l) && (a.m == b.m);
}

bool operator!=(const lm &a, const lm &b)
{
    return !(a == b);
}

double lerp(double x, double x0, double x1, double v0, double v1)
{
    return 1./(x1 - x0) * ((x1 - x)*v0 + (x - x0)*v1);
}
