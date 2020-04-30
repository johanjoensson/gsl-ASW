#include "augmented_fun.h"
#include "spherical_fun.h"
#include "schroedinger.h"

#include <numeric>
#include <fstream>

Augmented_function::Augmented_function(const int n_n, const lm l_n, const double kappa_n,
    const GSL::Vector& center_n, const Logarithmic_mesh& mesh_n)
 : n(n_n), l(l_n), kappa(kappa_n), radius(mesh_n.r_back()), center(center_n), mesh(mesh_n),
 val(mesh_n.size())
{}

double Augmented_function::operator()(const GSL::Vector& r) const
{
    double res = 0;
    GSL::Vector ri = r - center;
    size_t t = 1;
    if(ri.norm<double>() <= radius){
        while(mesh.r(t) < ri.norm<double>() && t < mesh.size()){
            t++;
        }
        // r[t - 1] < |ri| and r[t] > |ri|
        if(t < mesh.size()){
            res = lerp(ri.norm<double>(), mesh.r(t - 1), mesh.r(t), val[t - 1]*std::sqrt(mesh.drx(t-1))*mesh.drx(t-1), val[t]*std::sqrt(mesh.drx(t))*mesh.drx(t));
        }
    }
    return res;
}

double augmented_integral(const Augmented_function &a, const Augmented_function &b)
{
    Logarithmic_mesh mesh = a.mesh;
    // Make sure both funcions are centered on the same site
    if(a.center != b.center){
        return 0.;
    }

    std::vector<double> integrand(a.val.size(), 0);

    for(size_t i = 0; i < integrand.size(); i++){
    	integrand[i] = a.val[i]*b.val[i]*mesh.drx(i);
    }
//    return mesh.integrate(integrand);
    return mesh.integrate_simpson(integrand);
}

bool operator==(const Augmented_function &a, const Augmented_function &b)
{
    return a.n == b.n && a.l == b.l && a.center == b.center;
}

bool operator!=(const Augmented_function &a, const Augmented_function &b)
{
    return !(a == b);
}

Augmented_Hankel::Augmented_Hankel(const int n_n, const lm l_n, const double kappa_n,
    const GSL::Vector& center_n, const Logarithmic_mesh& mesh_n)
 : Augmented_function(n_n, l_n, kappa_n, center_n, mesh_n), EH()
{}

void Augmented_Hankel::update(std::vector<double>& v, const double en
    , const bool core)
{
    EH = en;
    size_t nodes = static_cast<size_t>(n - l.l - 1);
    size_t last = mesh.size() - 1, lastbutone = mesh.size() - 2;

    Hankel_function H(l);

    std::vector<double> l_init = {
        0,
        GSL::pow_int(mesh.r(1), l.l + 1)/std::sqrt(mesh.drx(1)),
        GSL::pow_int(mesh.r(2), l.l + 1)/std::sqrt(mesh.drx(2))};
    std::vector<double> r_init;
    if(core){
        r_init = {1e-8, 0.};
    }else{
        r_init = {GSL::pow_int(kappa, l.l + 1)*mesh.r(lastbutone)*H(kappa*mesh.r(lastbutone))/std::sqrt(mesh.drx(lastbutone)),
                  GSL::pow_int(kappa, l.l + 1)*mesh.r(last)*H(kappa*mesh.r(last))/std::sqrt(mesh.drx(last))};
    }
    Radial_Schroedinger_Equation_Central_Potential se(v, static_cast<size_t>(l.l), l_init, r_init, mesh, 1e-10);
    se.solve(nodes, EH);
    if(core){
        se.normalize();
    }else{
        double scale = 1;
        if(std::abs(se.psi().back()) > 1e-12){
            scale = r_init.back()/se.psi().back();
        }
        for(size_t i = 0; i < mesh.size(); i++){
            se.psi()[i] *= scale;
        }
    }
    val = se.psi();
    EH = se.e();
}

Augmented_Bessel::Augmented_Bessel(const int n_n, const lm l_n, const double kappa_n,
    const GSL::Vector& center_n, const Logarithmic_mesh& mesh_n)
 : Augmented_function(n_n, l_n, kappa_n, center_n, mesh_n), EJ()
{}

void Augmented_Bessel::update(std::vector<double>& v, const double en
    , const bool core)
{
    EJ = en;
    size_t nodes =  static_cast<size_t>(n - l.l - 1);
    int sign = 1;
    if(nodes % 2 == 1){
	    sign = -1;
    }
    size_t last = mesh.size() - 1, lastbutone = mesh.size() - 2;
    Bessel_function I(l);

    if(!core){
        std::vector<double> l_init = {
            0,
            GSL::pow_int(mesh.r(1), l.l+1)/std::sqrt(mesh.drx(1)),
            GSL::pow_int(mesh.r(2), l.l+1)/std::sqrt(mesh.drx(2))};

        std::vector<double> r_init = {
            sign*GSL::pow_int(kappa, -l.l)* I(kappa*mesh.r(lastbutone))
		        *mesh.r(lastbutone)/std::sqrt(mesh.drx(lastbutone)),
		    sign*GSL::pow_int(kappa, -l.l)*I(kappa*mesh.r(last))
		        *mesh.r(last)/std::sqrt(mesh.drx(last))};

        Radial_Schroedinger_Equation_Central_Potential se(v, static_cast<size_t>(l.l), l_init, r_init, mesh, 1e-10);
        se.solve(nodes, EJ);
        double scale = r_init.back()/se.psi().back();
        for(size_t i = 0; i < mesh.size(); i++){
            se.psi()[i] *= scale;
        }
        val = se.psi();
        EJ = se.e();
    }else{
        EJ = 0.;
        val = std::vector<double>(val.size(), 0.);
    }
}
