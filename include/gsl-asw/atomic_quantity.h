#ifndef ATOMIC_QUANTITY_H
#define ATOMIC_QUANTITY_H
#include <vector>
#include <functional>
#include <GSLpp/block.h>
#include <GSLpp/matrix.h>
#include "xc_func.h"
#include <numerical-mesh/numerical-mesh.h>

class Atomic_quantity;
class Potential;
class Density;

class Atomic_quantity{
protected:
    std::vector<size_t> mesh_lengths_m, lmax_m;
    bool spinpol_m;

public:
    Atomic_quantity() = default;
    Atomic_quantity(const Atomic_quantity&) = default;
    Atomic_quantity(Atomic_quantity&&) = default;
    virtual ~Atomic_quantity(){}

    Atomic_quantity(const std::vector<size_t>& mesh_lengths, const std::vector<size_t>& lmax, bool spinpol)
     : mesh_lengths_m(mesh_lengths), lmax_m(lmax), spinpol_m(spinpol)
    {}

    Atomic_quantity& operator=(const Atomic_quantity&) = default;
    Atomic_quantity& operator=(Atomic_quantity&&) = default;
    
    size_t n_at() const {return lmax_m.size();}
    
    const std::vector<size_t>& mesh_lengths() const { return mesh_lengths_m;}
    const std::vector<size_t>& lmax() const { return lmax_m;}
    
    bool spinpol() const {return spinpol_m;}
};

class Potential : public Atomic_quantity{

    GSL::Block atomic_m, hartree_m, exchange_correlation_pot_m, exchange_correlation_energy_m;
    Xc_func xc_fun_m;

    // double Xi0(const Site_t<3>& j, const double r);

    void initial_pot(const std::vector<Exponential_mesh<1, double>>& at_meshes, const std::vector<size_t>& zs, const Density& rho);
    void calc_atomic(const std::vector<Exponential_mesh<1, double>>& at_meshes, const std::vector<size_t>& zs);
    void calc_Hartree(const std::vector<Exponential_mesh<1, double>>& at_meshes, const Density& rho);
    void calc_XC(const Density& rho);

public:
    Potential() = default;
    Potential(const Potential&) = default;
    Potential(Potential&&) = default;
    ~Potential(){}

    Potential(const std::vector<Exponential_mesh<1, double>>& at_meshes, const std::vector<size_t>& zs, const std::vector<size_t>& l_max,
             const Density& rho, Xc_func xcf, bool spinpol = false);
    Potential& operator=(const Potential&) = default;
    Potential& operator=(Potential&&) = default;

    double MT0(){return 0;}

    GSL::Vector::View atomic(const size_t i);
    GSL::Vector::Const_View atomic(const size_t i) const;

    GSL::Matrix::View Hartree(const size_t i);
    GSL::Matrix::Const_View Hartree(const size_t i) const;
    GSL::Vector::View Hartree(const size_t i, const size_t l);
    GSL::Vector::Const_View Hartree(const size_t i, const size_t l) const;

    GSL::Matrix::View Vxc(const size_t i);
    GSL::Matrix::Const_View Vxc(const size_t i) const;
    GSL::Vector::View Vxc(const size_t i, const size_t l);
    GSL::Vector::Const_View Vxc(const size_t i, const size_t l) const;

    const GSL::Matrix tot(const size_t i) const;
    const GSL::Vector tot(const size_t i, const size_t l) const;
/**/
    GSL::Matrix::View exc(const size_t i);
    GSL::Matrix::Const_View exc(const size_t i) const;
    GSL::Vector::View exc(const size_t i, const size_t l);
    GSL::Vector::Const_View exc(const size_t i, const size_t l) const;

    void calc_pot(const std::vector<Exponential_mesh<1, double>>& at_meshes, const Density& rho);

};

class Density : public Atomic_quantity{
    GSL::Block valence_m, core_m;
public:
    Density() = default;
    Density(const Density&) = default;
    Density(Density&&) = default;
    ~Density(){}

    Density(const std::vector<Exponential_mesh<1, double>>& at_meshes, const std::vector<size_t>& l_max,
            bool spinpol = false);

    Density(const std::vector<size_t>& mesh_lengths, const std::vector<size_t> l_max,
        bool spinpol = false);

    GSL::Matrix::View valence(size_t i);
    GSL::Matrix::Const_View valence(size_t i) const;
    GSL::Vector::View valence(size_t i, size_t l);
    GSL::Vector::Const_View valence(size_t i, size_t l) const;

    GSL::Vector::View core(size_t i);
    GSL::Vector::Const_View core(size_t i) const;

    GSL::Vector tot(size_t i) const;

    Density& operator=(const Density&) = default;
    Density& operator=(Density&&) = default;

};
#endif // ATOMIC_QUANTITY_H
