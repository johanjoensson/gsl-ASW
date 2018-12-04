#ifndef XC_FUNC_H
#define XC_FUNC_H

#include <xc.h>
#include <vector>
#include <memory>
#include <functional>

enum XC_FUN{
    LDA,
    LDA_NM,
    UNDEFINED
};

class Xc_func{
private:
    std::unique_ptr<xc_func_type, std::function<void(xc_func_type*)>> fun;
public:
    Xc_func();
    Xc_func(XC_FUN xcf);
    Xc_func(Xc_func &xcf);
    Xc_func(Xc_func &&xcf);

    ~Xc_func();

    Xc_func& operator=(const Xc_func &xcf);
    Xc_func& operator=(Xc_func &&xcf);

    void set_xc(XC_FUN xcf);

    std::vector<double> exc(std::vector<double> rho);
};

#endif // XC_FUNC_H
