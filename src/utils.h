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
    lm(int l_n, int m_n):n(0), l(l_n), m(m_n){}
    lm():n(0), l(0), m(0){}
    int n;
    int l;
	int m;

    bool operator==(const lm &a) const
    {
        return (this->n == a.n) && (this->l == a.l) && (this->m == a.m);
    }

    bool operator!=(const lm &a) const
    {
        return !(*this == a);
    }
};

inline std::ostream& operator<<(std::ostream& os, const lm& l)
{
    return os << "(" << l.n << ", " << l.l << ", " << l.m <<")";
}

struct lm_hash{
    size_t operator()(const lm& l)
    {
        size_t res = 0;
        res ^= std::hash<int>()(l.n) + 0x9e3779b9 + (res<< 6) + (res>> 2);
        res ^= std::hash<int>()(l.l) + 0x9e3779b9 + (res<< 6) + (res>> 2);
        res ^= std::hash<int>()(l.m) + 0x9e3779b9 + (res<< 6) + (res>> 2);
        return res;
    }
};

template<class X = double, class Y = double>
double lerp(X x, X x0, X x1, Y y0, Y y1)
{
    return 1/(x1 - x0) * ((x1 - x)*y0 + (x - x0)*y1);
}

double calc_eta();

double calc_Rmax(const double vol, const double kappa, const lm &l, const double tol);
double calc_Kmax(const double vol, const double kappa, const lm &l, const double tol);


#endif // UTILS_H
