#ifndef UTILS_H
#define UTILS_H

enum spin{
    UP,
    DOWN,
    NON_COLLINEAR
};

struct lm {
	int l;
	int m;
};

double lerp(double x, double x0, double x1, double v0, double v1);

bool operator==(const lm &a, const lm &b);
bool operator!=(const lm &a, const lm &b);

#endif // UTILS_H
