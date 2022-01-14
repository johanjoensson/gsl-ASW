#ifndef SIMULATION_H
#define SIMULATION_H
#include <vector>
#include <map>
#include <functional>
#include "atomic_quantity.h"
#include "augmented_spherical_wave.h"
#include "crystal.h"
#include "site.h"
#include "xc_func.h"
#include "k-mesh.h"
#include "log_mesh.h"
#include "structure_const.h"
#include <numerical-mesh/numerical-mesh.h>
#include <GSLpp/matrix.h>
#include <GSLpp/matrix_complex.h>
#include <GSLpp/vector.h>


struct GSLVecCompare{
    bool operator()(const GSL::Vector& lhs, const GSL::Vector& rhs) const
    {
        return (lhs.norm() < rhs.norm());
    }
};

class Simulation{
    std::vector<double> kappas;
    std::vector<Hankel_container> Hs_m;
    std::vector<Bessel_container> Bs_m;
    Crystal_t<3, Atom> cryst;
    std::vector<Exponential_mesh<1, double>> at_meshes;
    Potential pot;
    Density n_m;
    std::vector<Augmented_spherical_wave> basis_valence;
    std::vector<Augmented_spherical_wave> basis_core;
    std::map<GSL::Vector::Const_View, std::pair<GSL::Matrix_Complex, GSL::Vector>, GSLVecCompare> k_eigenvals;
    std::vector<GSL::Matrix> XH1, XS1, XH2, XS2, XH3, XS3;
    Bloch_summed_structure_constant::Container B_m;

    void set_up_crystal();
    void set_up_basis();
    void set_up_augmented_functions();
    void set_up_potential(int xc_func_id);
    void init_augmented_functions();


    void add_states(const Site_t<3>& center, const double kappa);

    double X_H1(const Augmented_Hankel& Ht1, const Augmented_Hankel& Ht2, const Site_t<3>& at);
    double X_H2(const Augmented_Hankel& Ht1, const Augmented_Bessel& Jt2, const Site_t<3>& at);
    double X_H3(const Augmented_Bessel& Jt1, const Augmented_Bessel& Jt2, const Site_t<3>& at);
    double X_S1(const Augmented_Hankel& Ht1, const Augmented_Hankel& Ht2, const Site_t<3>& at);
    double X_S2(const Augmented_Hankel& Ht1, const Augmented_Bessel& Jt2, const Site_t<3>& at);
    double X_S3(const Augmented_Bessel& Jt1, const Augmented_Bessel& Jt2, const Site_t<3>& at);


    GSL::Complex H_element(const size_t i1, const size_t i2, const GSL::Vector& kp);
    GSL::Complex S_element(const size_t i1, const size_t i2, const GSL::Vector& kp);
public:
    Simulation(Crystal_t<3, Atom>&& crystal, int xc_func_id, const std::vector<double> kappas_n,
    std::function<double(const size_t, const double)> at_pot =
        [](const size_t Z, const double r)
        {
            return -2.*static_cast<double>(Z)/r;
        });

    void set_up_X_matrices();
    GSL::Matrix_Complex set_H(const GSL::Vector& kp);
    GSL::Matrix_Complex set_S(const GSL::Vector& kp);
    // const std::pair<GSL::Matrix_cx, GSL::Vector> calc_eigen(const GSL::Vector&) const;
    void calc_eigen(const GSL::Vector&);

    void add_eigvals(const GSL::Vector& kp, const GSL::Vector& eigvals);
    void print_eigvals(const K_mesh& k_mesh);

    double canonical_band(const lm l, const double kappa, const spin s, const GSL::Vector& kp);

};
#endif // SIMULATION_H
