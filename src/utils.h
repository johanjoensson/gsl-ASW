#ifndef UTILS_H
#define UTILS_H
#include <functional>
#include <iostream>

enum spin{
    UP,
    DOWN,
    NON_COLLINEAR
};

struct lm {
    int l;
	int m;
};

std::ostream& operator<<(std::ostream& os, const lm& l);

double lerp(double x, double x0, double x1, double v0, double v1);

bool operator==(const lm &a, const lm &b);
bool operator!=(const lm &a, const lm &b);

#endif // UTILS_H
