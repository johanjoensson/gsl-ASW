#include "xc_func.h"
#include <iostream>

Xc_func::Xc_func()
 : fun(nullptr)
{}


void Xc_func::set_xc(XC_FUN xcf)
{
    int xc_func = 0;
    int xc_spin = 0;

    if(fun.get() == nullptr){
        fun = std::unique_ptr<xc_func_type>(new xc_func_type);
    }

    switch(xcf) {
        case LDA :
            xc_func = XC_LDA_X;
            xc_spin = XC_POLARIZED;
            break;
        case LDA_NM :
            xc_func = XC_LDA_X;
            xc_spin = XC_UNPOLARIZED;
            break;
        default :
            break;
    }
    if(xc_func_init(fun.get(), xc_func, xc_spin) != 0){
        std::cerr << "Functional " << xc_func << " not found." << std::endl;
    }
}

Xc_func::Xc_func(XC_FUN xcf)
 : Xc_func()
{
    if(xcf != UNDEFINED){
        fun = std::unique_ptr<xc_func_type>(new xc_func_type);
        this->set_xc(xcf);
    }
}

Xc_func::Xc_func(Xc_func &xcf)
 : fun(new xc_func_type)
{
    *this->fun = *xcf.fun;
}

Xc_func::Xc_func(Xc_func &&xcf)
 : fun(nullptr)
{
    std::swap(this->fun, xcf.fun);
}

Xc_func& Xc_func::operator=(const Xc_func &xcf)
{
    if(this->fun == nullptr){
        this->fun = std::unique_ptr<xc_func_type>(new xc_func_type);
    }
    *this->fun = *xcf.fun;

    return *this;
}


Xc_func& Xc_func::operator=(Xc_func &&xcf)
{
    std::swap(this->fun, xcf.fun);

    return *this;
}

Xc_func::~Xc_func()
{
    if(fun.get() != nullptr){
        xc_func_end(fun.get());
    }
}

std::vector<double> Xc_func::exc(std::vector<double> rho)
{
    std::vector<double> res(rho.size(), 0.);

    switch(fun->info->family) {
        case XC_FAMILY_LDA:
            xc_lda_exc(fun.get(), rho.size(), &(*rho.begin()), &(*res.begin()));
            break;
        default:
            break;
    }


    return res;
}
