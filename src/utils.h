#ifndef UTILS_H
#define UTILS_H
/***************************************************************************//**
* Data type for representing spin.\n
*******************************************************************************/
enum spin{
    UP,
    DOWN,
    NON_COLLINEAR
};

/***************************************************************************//**
* Data type for representing combined orbital and magnetic quantum numbers.\n
*******************************************************************************/
struct lm {
    int l;
	int m;
};

//! Linear interpolation between (x0,v0) and (x1,v1) for x0 <= x <= x1
double lerp(double x, double x0, double x1, double v0, double v1);

bool operator==(const lm &a, const lm &b);
bool operator!=(const lm &a, const lm &b);

#endif // UTILS_H
