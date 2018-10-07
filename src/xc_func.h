#ifndef XC_FUNC_H
#define XC_FUNC_H

#include <xc.h>
#include <vector>

enum XC_FUN{
    LDA,
    LDA_NM,
    UNDEFINED
};

class Xc_func{
private:
    xc_func_type *fun;
public:
    Xc_func();
    Xc_func(XC_FUN xcf);

    ~Xc_func();
    void set_xc(XC_FUN xcf);

    std::vector<double> exc(std::vector<double> rho);
};

#endif // XC_FUNC_H
