#include <xc_func.h>
#include <iostream>
#include <algorithm>


Xc_func::Xc_func(int xcf_id)
 : Xc_func()
{
    if(xc_func_init(&fun, xcf_id, XC_UNPOLARIZED) != 0){
        fprintf(stderr, "Functional '%d' not found\n", xcf_id);
    }
}

GSL::Vector Xc_func::exc(GSL::Vector::Const_View rho)
{
    GSL::Vector Exc(rho.size());

    switch(fun.info->family) {
        case XC_FAMILY_LDA:
            xc_lda_exc(&fun, static_cast<int>(rho.size()), rho.gsl_data()->data, Exc.gsl_data()->data);
            break;
        default:
            break;
    }
    /* LibXC uses Hartree atomic units, we use Rydberg. 1 Ha = 2 Ry. */
    return 2*Exc;
}

GSL::Vector Xc_func::vxc(GSL::Vector::Const_View rho)
{
    GSL::Vector Vxc(rho.size());

    switch(fun.info->family) {
        case XC_FAMILY_LDA:
            xc_lda_vxc(&fun, static_cast<int>(rho.size()), rho.gsl_data()->data, Vxc.gsl_data()->data);
            break;
        default:
            break;
    }
    /* LibXC uses Hartree atomic units, we use Rydberg. 1 Ha = 2 Ry. */
    return 2*Vxc;
}

std::pair<GSL::Vector, GSL::Vector> Xc_func::exc_vxc(GSL::Vector::Const_View rho)
{
    GSL::Vector Vxc(rho.size()), Exc(rho.size());

    switch(fun.info->family) {
        case XC_FAMILY_LDA:
            xc_lda_exc_vxc(&fun, static_cast<int>(rho.size()), rho.gsl_data()->data, Exc.gsl_data()->data, Vxc.gsl_data()->data);
            break;
        default:
            break;
    }
    /* LibXC uses Hartree atomic units, we use Rydberg. 1 Ha = 2 Ry. */
    return {2*Exc, 2*Vxc};
}
