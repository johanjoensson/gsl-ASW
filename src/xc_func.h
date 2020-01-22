#ifndef XC_FUNC_H
#define XC_FUNC_H

#include <xc.h>
#include <vector>
#include <memory>
#include <functional>

enum XC_FUN{
    LDA,
    LDA_NM,
    NONE,
    UNDEFINED
};

class Xc_func{
private:
    std::unique_ptr<xc_func_type, std::function<void(xc_func_type*)>> fun;
public:
    Xc_func() : fun() {};
    Xc_func(const Xc_func&) = default;    
    Xc_func(Xc_func&& other) = default;
    Xc_func(XC_FUN xcf);
    ~Xc_func() = default;

    void set_xc(XC_FUN xcf);

    Xc_func& operator=(const Xc_func& other)
    {
        this->fun.reset(new xc_func_type (*other.fun));
        return *this;
    }
    Xc_func& operator=(Xc_func&& other)
    {
        this->fun = std::move(other.fun);
        return *this;
    }

    std::vector<double> exc(std::vector<double> rho);
};

#endif // XC_FUNC_H
