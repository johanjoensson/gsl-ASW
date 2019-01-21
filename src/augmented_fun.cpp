#include "augmented_fun.h"
#include "spherical_fun.h"
#include "numerov_solver.h"

#include <numeric>
#include <fstream>

Augmented_function::Augmented_function()
 : n(), l(), kappa(), radius(), center(), mesh(), val()
{}

Augmented_function::Augmented_function(const Augmented_function &a)
 : n(a.n), l(a.l), kappa(a.kappa), radius(a.radius), center(a.center), mesh(a.mesh), val(a.val)
{}

Augmented_function::Augmented_function(Augmented_function &&a)
 : Augmented_function()
{
    std::swap(n, a.n);
    std::swap(l, a.l);
    std::swap(kappa, a.kappa);
    std::swap(radius, a.radius);

    std::swap(center, a.center);
    std::swap(mesh, a.mesh);
    std::swap(val, a.val);
}

Augmented_function::Augmented_function(const int n_n, const lm l_n, const double kappa_n,
    const GSL::Vector& center_n, const Logarithmic_mesh& mesh_n)
 : n(n_n), l(l_n), kappa(kappa_n), radius(mesh_n.r.back()), center(center_n), mesh(mesh_n),
 val(mesh_n.r.size())
{}

Augmented_function::~Augmented_function()
{}

double Augmented_function::operator()(const GSL::Vector& r)
{
    double res = 0;
    GSL::Vector ri = r - center;
    size_t t = 1;
    if(ri.norm() <= radius){
        while(mesh.r[t] < ri.norm() && t < mesh.r.size()){
            t++;
        }
        // r[t - 1] < |ri| and r[t] > |ri|
        if(t < mesh.r.size()){
            res = lerp(ri.norm(), mesh.r[t - 1], mesh.r[t], val[t - 1]*mesh.drx[t-1], val[t]*mesh.drx[t]);
        }
    }
    return res;
}

Augmented_function& Augmented_function::operator=(const Augmented_function& a)
{
    n = a.n;
    l = a.l;
    kappa = a.kappa;
    radius = a.radius;
    center = a.center;
    mesh = a.mesh;
    val = a.val;

    return *this;

}

Augmented_function& Augmented_function::operator=(Augmented_function&& a)
{
    std::swap(n, a.n);
    std::swap(l, a.l);
    std::swap(kappa, a.kappa);
    std::swap(radius, a.radius);
    std::swap(center, a.center);
    std::swap(mesh, a.mesh);
    std::swap(val, a.val);

    return *this;

}

double augmented_integral(const Augmented_function &a, const Augmented_function &b)
{
    Logarithmic_mesh mesh = a.mesh;
    // Make sure both funcions are centered on the same site
    if(a.center != b.center){
        return 0.;
    }

//    res = 1./3 * res;
    std::vector<double> integrand(a.val.size(), 0);
    // integrand = r*f(r)*r*g(r) = r^2*f(r)*g(r)
    for(size_t i = 0; i < integrand.size(); i++){
    	integrand[i] = a.val[i]*b.val[i];
    }
    return mesh.integrate(integrand);
}

bool operator==(const Augmented_function &a, const Augmented_function &b)
{
    return a.n == b.n && a.l == b.l && a.center == b.center;
}

bool operator!=(const Augmented_function &a, const Augmented_function &b)
{
    return !(a == b);
}

Augmented_Hankel::Augmented_Hankel()
 : Augmented_function(), EH()
{}

Augmented_Hankel::Augmented_Hankel(const Augmented_Hankel& a)
 : Augmented_function(a), EH(a.EH)
{}

Augmented_Hankel::Augmented_Hankel(Augmented_Hankel&& a)
 : Augmented_function(), EH()
{
    std::swap(n, a.n);
    std::swap(l, a.l);
    std::swap(kappa, a.kappa);
    std::swap(radius, a.radius);

    std::swap(center, a.center);
    std::swap(mesh, a.mesh);
    std::swap(val, a.val);

    std::swap(EH, a.EH);
}

Augmented_Hankel::Augmented_Hankel(const int n_n, const lm l_n, const double kappa_n,
    const GSL::Vector& center_n, const Logarithmic_mesh& mesh_n)
 : Augmented_function(n_n, l_n, kappa_n, center_n, mesh_n), EH()
{}

Augmented_Hankel::~Augmented_Hankel()
{}

void Augmented_Hankel::update(std::vector<double>& v, const double en
    , const bool core)
{
    EH = en;
    Numerov_solver sol;
    int nodes = n - l.l - 1;
    int sign = 1;
    if(nodes % 2 == 1){
	    sign = -1;
    }
    size_t last = mesh.r.size() - 1, lastbutone = mesh.r.size() - 2;

    Hankel_function H(l);

    std::vector<double> l_init = {
        0,
        GSL::pow_int(mesh.r[1], l.l + 1),
        GSL::pow_int(mesh.r[2], l.l + 1)};
    std::vector<double> r_init;
    if(core){
        r_init = {0., 0.};
    }else{
        r_init = {sign*GSL::pow_int(kappa, l.l + 1)*mesh.r[lastbutone]*H(kappa*mesh.r[lastbutone]),
                  sign*GSL::pow_int(kappa, l.l + 1)*mesh.r[last]*H(kappa*mesh.r[last])};
    }
    val = sol.solve(mesh, v, l_init, r_init, EH, nodes);
}

Augmented_Hankel& Augmented_Hankel::operator=(const Augmented_Hankel& a)
{
    n = a.n;
    l = a.l;
    kappa = a.kappa;
    radius = a.radius;

    center = a.center;
    mesh = a.mesh;
    val = a.val;

    EH = a.EH;

    return *this;
}

Augmented_Hankel& Augmented_Hankel::operator=(Augmented_Hankel&& a)
{
    std::swap(n, a.n);
    std::swap(l, a.l);
    std::swap(kappa, a.kappa);
    std::swap(radius, a.radius);

    std::swap(center, a.center);
    std::swap(mesh, a.mesh);
    std::swap(val, a.val);

    std::swap(EH, a.EH);

    return *this;
}

Augmented_Bessel::Augmented_Bessel()
 : Augmented_function(), EJ()
{}

Augmented_Bessel::Augmented_Bessel(const Augmented_Bessel& a)
 : Augmented_function(a), EJ(a.EJ)
{
    EJ = a.EJ;
}

Augmented_Bessel::Augmented_Bessel(Augmented_Bessel&& a)
 : Augmented_function(), EJ()
{
    std::swap(n, a.n);
    std::swap(l, a.l);
    std::swap(kappa, a.kappa);
    std::swap(radius, a.radius);

    std::swap(center, a.center);
    std::swap(mesh, a.mesh);
    std::swap(val, a.val);

    std::swap(EJ, a.EJ);
}

Augmented_Bessel::Augmented_Bessel(const int n_n, const lm l_n, const double kappa_n,
    const GSL::Vector& center_n, const Logarithmic_mesh& mesh_n)
 : Augmented_function(n_n, l_n, kappa_n, center_n, mesh_n), EJ()
{}

Augmented_Bessel::~Augmented_Bessel()
{}

void Augmented_Bessel::update(std::vector<double>& v, const double en
    , const bool core)
{
    EJ = en;
    Numerov_solver sol;
    int nodes =  n - l.l - 1;
    int sign = 1;
    if(nodes % 2 == 1){
	    sign = -1;
    }
    size_t last = mesh.r.size() - 1, lastbutone = mesh.r.size() - 2;
    Bessel_function i(l);

    if(!core && nodes >= 0){
        std::vector<double> l_init = {
            0,
            GSL::pow_int(mesh.r[1], l.l+1),
            GSL::pow_int(mesh.r[2], l.l+1)};

        std::vector<double> r_init = {
                sign*GSL::pow_int(kappa, -l.l)*i(kappa*mesh.r[lastbutone])
		        *mesh.r[lastbutone],
		sign*GSL::pow_int(kappa, -l.l)*i(kappa*mesh.r[last])
		        *mesh.r[last]};

        val = sol.solve(mesh, v, l_init, r_init, EJ, nodes);
    }else{
        EJ = 0.;
        val = std::vector<double>(val.size(), 0.);
    }
}

Augmented_Bessel& Augmented_Bessel::operator=(const Augmented_Bessel& a)
{
    n = a.n;
    l = a.l;
    kappa = a.kappa;
    radius = a.radius;
    EJ = a.EJ;

    center = a.center;
    mesh = a.mesh;
    val = a.val;

    return *this;
}

Augmented_Bessel& Augmented_Bessel::operator=(Augmented_Bessel&& a)
{
    std::swap(n, a.n);
    std::swap(l, a.l);
    std::swap(kappa, a.kappa);
    std::swap(radius, a.radius);

    std::swap(center, a.center);
    std::swap(mesh, a.mesh);
    std::swap(val, a.val);

    std::swap(EJ, a.EJ);

    return *this;
}
