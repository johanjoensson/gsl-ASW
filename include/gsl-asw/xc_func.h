#ifndef XC_FUNC_H
#define XC_FUNC_H
#include <xc.h>
#include <vector>
#include <memory>
#include <functional>
#include <utility>
#include <GSLpp/vector.h>

class Xc_func{
private:
    xc_func_type fun;

public:
    Xc_func() = default;
    Xc_func(const Xc_func&) = default;
    Xc_func(Xc_func&& other) = default;
    Xc_func(int xcf_id);
    ~Xc_func() {xc_func_end(&fun);}

    Xc_func& operator=(const Xc_func& other) = default;
    Xc_func& operator=(Xc_func&& other) = default;

    GSL::Vector exc(GSL::Vector::Const_View rho);
    GSL::Vector vxc(GSL::Vector::Const_View rho);
    std::pair<GSL::Vector, GSL::Vector> exc_vxc(GSL::Vector::Const_View rho);

};

#endif // XC_FUNC_H
