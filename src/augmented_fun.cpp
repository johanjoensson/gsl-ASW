#include "augmented_fun.h"
#include "spherical_fun.h"
#include "numerov_solver.h"

Augmented_function::Augmented_function()
 : n(), l(), kappa(), radius(), center(), mesh(), val()
{}

Augmented_function::Augmented_function(const Augmented_function &a)
{
    n = a.n;
    l = a.l;
    kappa = a.kappa;
    radius = a.radius;
    center = a.center;
    mesh = a.mesh;
    val = a.val;
}

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

Augmented_function::Augmented_function(const int n, const lm l, const double kappa,
    const GSL::Vector center, const Logarithmic_mesh mesh)
 : n(n), l(l), kappa(kappa), radius(mesh.r.back()), center(center), mesh(mesh),
 val(mesh.r.size())
{
}

Augmented_function::~Augmented_function()
= default;

double Augmented_function::operator()(const GSL::Vector r)
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
            res = lerp(ri.norm(), mesh.r[t - 1], mesh.r[t], val[t - 1], val[t]);
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

bool operator==(const Augmented_function &a, const Augmented_function &b)
{
    return a.l == b.l && a.center == b.center;
}

bool operator!=(const Augmented_function &a, const Augmented_function &b)
{
    return !(a == b);
}

Augmented_Hankel::Augmented_Hankel()
 : Augmented_function(), EH()
{}

Augmented_Hankel::Augmented_Hankel(const Augmented_Hankel& a)
 : Augmented_function(a)
{
    n = a.n;
    l = a.l;
    kappa = a.kappa;
    radius = a.radius;
    center = a.center;
    mesh = a.mesh;
    val = a.val;

    EH = a.EH;

}

Augmented_Hankel::Augmented_Hankel(Augmented_Hankel&& a)
 : Augmented_function(a), EH()
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

Augmented_Hankel::Augmented_Hankel(const int n, const lm l, const double kappa,
    const GSL::Vector center, const Logarithmic_mesh mesh)
 : Augmented_function(n, l, kappa, center, mesh), EH()
{}

Augmented_Hankel::~Augmented_Hankel()
{}

void Augmented_Hankel::update(std::vector<double> v, const double en
    , const bool core)
{
    EH = en;
    Numerov_solver sol;
    int nodes = n - l.l - 1;

    std::vector<double> l_init = {GSL::pow_int(mesh.r[0], l.l+1),
        GSL::pow_int(mesh.r[1], l.l+1),
        GSL::pow_int(mesh.r[2], l.l+1)};
    std::vector<double> r_init;
    if(core){
        r_init = {1e-16, 0.};
    }else{
        r_init = {GSL::pow_int(kappa, l.l + 1)*
            real_spherical_hankel(l, kappa*mesh.r[mesh.r.size()-1]).val,
                  GSL::pow_int(kappa, l.l+1)*
            real_spherical_hankel(l, kappa*mesh.r[mesh.r.size()-2]).val};
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
 : Augmented_function(a)
{
    n = a.n;
    l = a.l;
    kappa = a.kappa;
    radius = a.radius;
    EJ = a.EJ;

    center = a.center;
    mesh = a.mesh;
    val = a.val;
}

Augmented_Bessel::Augmented_Bessel(Augmented_Bessel&& a)
 : Augmented_function(a), EJ()
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

Augmented_Bessel::Augmented_Bessel(const int n, const lm l, const double kappa,
    const GSL::Vector center, const Logarithmic_mesh mesh)
 : Augmented_function(n, l, kappa, center, mesh), EJ()
{}

Augmented_Bessel::~Augmented_Bessel()
{}

void Augmented_Bessel::update(std::vector<double> v, const double en
    , const bool core)
{
    EJ = en;
    Numerov_solver sol;
    int nodes = n - l.l - 1;

    if(!core){
        std::vector<double> l_init = {GSL::pow_int(mesh.r[0], l.l+1),
            GSL::pow_int(mesh.r[1], l.l+1),
            GSL::pow_int(mesh.r[2], l.l+1)};

        std::vector<double> r_init = {GSL::pow_int(1./kappa, l.l)*
            GSL::bessel_jn(l.l,kappa*mesh.r[mesh.r.size() - 1]).val,
                  GSL::pow_int(1./kappa, l.l)*
            GSL::bessel_jn(l.l, kappa*mesh.r[mesh.r.size() - 2]).val};

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
