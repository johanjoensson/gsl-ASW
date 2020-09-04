#ifndef UTILS_H
#define UTILS_H
#include <functional>
#include <iostream>

enum spin{
    UP,
    DOWN,
    NON_COLLINEAR
};

/*******************************************************************************
 Composite quantum numbers, (n, l, m)
*******************************************************************************/
struct lm {
    lm(const int& l_n, const int& m_n):n(0), l(l_n), m(m_n){}
    lm(const int& n_n, const int& l_n, const int& m_n):n(n_n), l(l_n), m(m_n){}
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

    /***************************************************************************
    Ordering according to the Aufbauprinzip i.e., first according to n + l,
    if n + l equal, order according to n. If n are equal order according to
    m (ordering according to m is not part of the Aufbauprinzip).
    E.g. [1s, 2s, 2p, 3s, 3p, 4s, 3d, 4p, ...]
       /
    1s
       /   /
    2s  2p
       /   /   /
    3s  3p  3d
       /   /
    4s  4p  4d  4f
    .
    .
    .
    ***************************************************************************/
    bool operator<(const lm& a) const
    {
        if(this->n + this->l == a.n + a.l){
            if(this->n == a.n){
                return (this->m < a.m);
            }else{
                return (this->n < a.n);
            }
        }else{
            return (this->n + this->l < a.n + a.l);
        }
        return false;
    }

    bool operator>(const lm& a) const
    {
        if(this->n + this->l == a.n + a.l){
            if(this->n == a.n){
                return (this->m > a.m);
            }else{
                return (this->n > a.n);
            }
        }else{
            return (this->n + this->l > a.n + a.l);
        }
        return false;
    }

    bool operator>=(const lm& a) const
    {
        return !(*this < a);
    }

    bool operator<=(const lm& a) const
    {
        return !(*this > a);
    }

    lm& operator++()
    {
        if(++this->m > this->l){
            if(++this->l > this->n - 1){
                ++this->n;
                this->l = 0;
            }
            this->m = -this->l;
        }
        return *this;
    }

    lm& operator--()
    {
        if(--this->m < -this->l){
            if(--this->l < 0){
                --this->n;
                this->l = this->n - 1;
            }
            this->m = this->l;
        }
        return *this;
    }

    lm operator++(int)
    {
        lm res = *this;
        ++*this;
        return res;
    }

    lm operator--(int)
    {
        lm res = *this;
        --*this;
        return res;
    }

    lm& operator+=(int a)
    {
        int tmp = a;
        while(tmp > this->l - this->m){
            tmp -= this->l - this->m + 1;
            this->l += 1;
            if(this->l >= this->n){
                this->l = 0;
                this->n++;
            }
            this->m = -this->l;
        }
        this->m += tmp;
        return *this;
    }

    lm& operator-=(int a)
    {
        int tmp = a;
        while(tmp > this->m + this->l){
            tmp -= this->l + this->m + 1;
            this->l -= 1;
            if(this->l < 0){
                this->n--;
                this->l = this->n - 1;
            }
            this->m = this->l;
        }
        this->m -= tmp;
        return *this;
    }

    lm operator+(int a) const
    {
        lm res = *this;
        return (res += a);
    }

    lm operator-(int a) const
    {
        lm res = *this;
        return (res -= a);
    }
    std::string to_string() const
    {
        std::string res = "";
        res += "(" + std::to_string(this->n);
        res += ", " + std::to_string(this->l);
        res += ", " + std::to_string(this->m) + ")";
        return res;
    }
};

inline std::ostream& operator<<(std::ostream& os, const lm& l)
{
    return os << l.to_string();
}

template<class X = double, class Y = double>
double lerp(X x, X x0, X x1, Y y0, Y y1)
{
    return 1/(x1 - x0) * ((x1 - x)*y0 + (x - x0)*y1);
}

double calc_eta(const double volume);

double calc_Rmax(const double vol, const double kappa, const lm &l, const double tol);
double calc_Kmax(const double vol, const double kappa, const lm &l, const double tol);


namespace std {
	template<>
    struct hash<lm>{
        size_t operator()(const lm& l)
        {
            size_t res = 0;
            res ^= std::hash<int>()(l.n) + 0x9e3779b9 + (res<< 6) + (res>> 2);
            res ^= std::hash<int>()(l.l) + 0x9e3779b9 + (res<< 6) + (res>> 2);
            res ^= std::hash<int>()(l.m) + 0x9e3779b9 + (res<< 6) + (res>> 2);
            return res;
        }
    };
}
#endif // UTILS_H
