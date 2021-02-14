#ifndef ASW_BASIS_H
#define ASW_BASIS_H

#include <GSLpp/matrix.h>
#include <GSLpp/vector.h>

class ASW_basis {
public:

    ASW_basis() = delete;
    ASW_basis(const ASW_basis&) = default;
    ASW_basis(ASW_basis&&) = default;

    ASW_basis(const Crystal<3, Atom>& c) : cryst_m{c} {}

    ~ASW_basis() = default;

    ASW_basis& operator=(const ASW_basis&) = default;
    ASW_basis& operator=(ASW_basis&&) = default;

protected:
private:
    const Crystal_t<3, Atom>& cryst_m;
    std::vector<double> linearization_energies_m;
    std::vector<Hankel_container> Hankel_m;
    std::vector<Bessel_container> Bessel_m;
    std::vector<GSL::Matrix> XH1, XS1, XH2, XS2, XH3, XS3;
    Bloch_summed_structure_constant::Container B_m;

    double X_H1(const Augmented_Hankel& Ht1, const Augmented_Hankel& Ht2, const Site_t<3>& at);
    double X_H2(const Augmented_Hankel& Ht1, const Augmented_Bessel& Jt2, const Site_t<3>& at);
    double X_H3(const Augmented_Bessel& Jt1, const Augmented_Bessel& Jt2, const Site_t<3>& at);
    double X_S1(const Augmented_Hankel& Ht1, const Augmented_Hankel& Ht2, const Site_t<3>& at);
    double X_S2(const Augmented_Hankel& Ht1, const Augmented_Bessel& Jt2, const Site_t<3>& at);
    double X_S3(const Augmented_Bessel& Jt1, const Augmented_Bessel& Jt2, const Site_t<3>& at);

    GSL::Matrix_cx hamiltonian(GSL::Vector kp) const;
    GSL::Matrix_cx overlap(GSL::Vector kp) const;
};

#endif // ASW_BASIS_H
