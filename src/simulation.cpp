#include "simulation.h"
#include <iostream>
#include <stdexcept>
#include <set>
#include <algorithm>
#include "atom.h"
#include "utils.h"
#include "envelope_fun.h"
#include "structure_const.h"
#include "GSLpp/special_functions.h"
#include "GSLpp/eigen.h"
#include "GSLpp/linalg.h"

Simulation::Simulation()
 : cryst(), pot(), n(), basis_valence(), basis_core(), H(), S(), XH1(), XS1(),
   XH2(), XS2(), XH3(), XS3()
{}

Simulation::Simulation(const Crystal& crystal, const XC_FUN func, const double kappa)
 : cryst(crystal), pot(), n(), basis_valence(), basis_core(), H(), S(), XH1(),
   XS1(), XH2(), XS2(), XH3(), XS3()
{
    // Calculate MT radii
    double at_vol = 0;
    size_t l_low = 2;
    std::vector<std::vector<Atom>> nn = cryst.calc_nearest_neighbours();
    for(size_t i = 0; i < cryst.atoms.size(); i++){
        cryst.atoms[i].set_MT(nn[i][0].get_pos().norm()/2);
		at_vol += GSL::pow_int(cryst.atoms[i].get_MT(), 3);
    }
	at_vol *= 4*M_PI/3;
    // Set AS radii, atomic sphere volumes should add up to the crystal volume
	for(size_t i = 0; i < cryst.atoms.size(); i++){
		cryst.atoms[i].set_AS( std::cbrt(cryst.volume/at_vol) *
		cryst.atoms[i].get_MT());
	}
    // Set logarithmic meshes for the atoms and count number of electrons
    size_t nel = 0;
	for(size_t i = 0; i < cryst.atoms.size(); i++){
        cryst.atoms[i].mesh = Logarithmic_mesh(cryst.atoms[i].get_AS(),
            500 + static_cast<size_t>(std::ceil(
                1000*(static_cast<double>(cryst.atoms[i].get_Z())/118.)
            ))
            );
        nel += cryst.atoms[i].get_Z();
    }
    // Set up initial potential
    pot = Potential(cryst.atoms);
    pot.set_xc_fun(func);
    pot.initial_pot(nel, cryst.volume);
    std::vector<Augmented_spherical_wave> tmp;
    // Divide electrons into core and valence states
    for(Atom& at: cryst.atoms){
        add_states(at, kappa);
        if(at.get_Z() > 20){
            l_low = 3;
        }
    }
    std::cout << "# valence states " << basis_valence.size() << std::endl;

    if(basis_core.size() >= nel){
        throw std::runtime_error("Number of core states exceeds or equals number of electrons in the system!");
    }
    std::cout << "Setting up core states." << std::endl;
    for(Augmented_spherical_wave& wave : basis_core){
        wave.core_state = true;
        wave.set_up(pot);
    }
    std::cout << "Setting up valence states." << std::endl;
    for(Augmented_spherical_wave& wave: basis_valence){
        wave.core_state = false;
        wave.set_up(pot);
    }
    size_t N = basis_valence.size();
    H = GSL::Complex_Matrix(N, N);
    S = GSL::Complex_Matrix(N, N);
    // Reserve space for each atom ad each 0 <= l <= l_low, -l <= m <= l
    // XH1 = GSL::Matrix(cryst.atoms.size(), (l_low + 1)*(l_low + 1));
    // XS1 = GSL::Matrix(cryst.atoms.size(), (l_low + 1)*(l_low + 1));
    // XH2 = GSL::Matrix(cryst.atoms.size(), (l_low + 1)*(l_low + 1));
    // XS2 = GSL::Matrix(cryst.atoms.size(), (l_low + 1)*(l_low + 1));
    XH1 = GSL::Matrix(cryst.atoms.size(),(l_low + 1)*(l_low + 1));
    XS1 = GSL::Matrix(cryst.atoms.size(),(l_low + 1)*(l_low + 1));
    XH2 = GSL::Matrix(cryst.atoms.size(),(l_low + 1)*(l_low + 1));
    XS2 = GSL::Matrix(cryst.atoms.size(),(l_low + 1)*(l_low + 1));
    // Reserve space for each atom and each 0<= l <= l_int = l_low + 1, -l <= m <= l
    XH3 = GSL::Matrix(cryst.atoms.size(),(l_low + 2)*(l_low + 2));
    XS3 = GSL::Matrix(cryst.atoms.size(),(l_low + 2)*(l_low + 2));
}

struct nl_comp {
    bool operator() (const std::pair<int, int>& lhs, const std::pair<int, int>& rhs)
    {
        if(lhs.first + lhs.second == rhs.first + rhs.second){
            return lhs.first < rhs.first;
        }else{
            return (lhs.first + lhs.second) < (rhs.first + rhs.second);
        }
    }
};

void Simulation::add_states(const Atom& at, const double kappa)
{
    std::set<std::pair<int, int>, nl_comp> p({{1, 0}});
    int n_s = 1, ln = 0;
    size_t nel = 0;
    size_t Z = at.get_Z();
    while(nel < Z){
        for(int l = 0; l < n_s; l++){
            p.insert({n_s+1, l});
            nel += static_cast<size_t>(2*(2*l + 1));
        }
        p.insert({n_s+1, n_s});
        nel += static_cast<size_t>(2*(2*n_s + 1));
        n_s++;
    }
    nel = 0;
    for(std::pair<int, int> tmp : p){
        n_s = tmp.first;
        ln = tmp.second;
        if(nel + static_cast<size_t>(2*n_s*n_s) >= Z){
            std::cout << "valence : ( " << n_s << " " << ln << " )" << std::endl;
            for(int m = -ln; m <= ln; m++){
                basis_valence.push_back(Augmented_spherical_wave(kappa, n_s, lm {ln, m}, UP, at, cryst.atoms));
                basis_valence.push_back(Augmented_spherical_wave(kappa, n_s, lm {ln, m}, DOWN, at, cryst.atoms));
            }
        }else{
            std::cout << "core : ( " << n_s << " " << ln << " )" << std::endl;
            for(int m = -ln; m <= ln; m++){
                basis_core.push_back(Augmented_spherical_wave(kappa, n_s, lm {ln, m}, UP, at, cryst.atoms));
                basis_core.push_back(Augmented_spherical_wave(kappa, n_s, lm {ln, m}, DOWN, at, cryst.atoms));
            }
        }
        nel += static_cast<size_t>(2*(2*ln + 1));
    }
}

double X_H1(const Augmented_Hankel& Ht1, const Augmented_Hankel& Ht2, const Atom& at)
{
    double res = 0;
    const Envelope_Hankel H1(at, Ht1.l, Ht1.kappa);
    const Envelope_Hankel H2(at, Ht2.l, Ht2.kappa);
    res += Ht2.EH*augmented_integral(Ht1, Ht2);
    res += -Ht2.kappa*Ht2.kappa*off_atomic_integral(H1, H2);
    return res;
}

double X_H2(const Augmented_Hankel& Ht1, const Augmented_Bessel& Jt2, const Atom& at)
{
    double res = 0;
    const Envelope_Hankel H1(at, Ht1.l, Ht1.kappa);
    const Envelope_Bessel J2(at, Jt2.l, Jt2.kappa);
    res += Jt2.EJ*augmented_integral(Ht1, Jt2);
    res -= -Jt2.kappa*Jt2.kappa*atomic_integral(H1, J2);
    if(Ht1.kappa*Ht1.kappa != Jt2.kappa*Jt2.kappa){
        res += -Jt2.kappa*Jt2.kappa/(-Jt2.kappa*Jt2.kappa - -Ht1.kappa*Ht1.kappa);
    }
    return res;
}

double X_H3(const Augmented_Bessel& Jt1, const Augmented_Bessel& Jt2, const Atom& at)
{
    double res = 0;
    const Envelope_Bessel J1(at, Jt1.l, Jt1.kappa);
    const Envelope_Bessel J2(at, Jt2.l, Jt2.kappa);
    res += Jt2.EJ*augmented_integral(Jt1, Jt2);
    res -= -Jt2.kappa*Jt2.kappa*atomic_integral(J1, J2);
    return res;
}

GSL::Complex Simulation::H_element(const size_t i1, const size_t i2, const GSL::Vector& kp)
{
    const Augmented_spherical_wave w1 = basis_valence[i1];
    const Augmented_spherical_wave w2 = basis_valence[i2];
    Bloch_summed_structure_constant B1, B2;

    int l_low_1 = 2, l_low_2 = 2, l_low  = 2;
    if(w1.center.Z > 20){
        l_low_1 = 3;
    }
    if(w2.center.Z > 20){
        l_low_2 = 3;
    }
    GSL::Complex res(0,0);
    const GSL::Vector tau_ij(w1.center.pos - w2.center.pos);
    size_t at_i = 0, at_j = 0;
    for(size_t i = 0; i < cryst.atoms.size(); i++){
        if(cryst.atoms[i] == w1.center){
            at_i = i;
        }
        if(cryst.atoms[i] == w2.center){
            at_j = i;
        }
    }
    size_t l_i = static_cast<size_t>(w1.l.l*(w1.l.l + 1) + w1.l.m);
    size_t l_j = static_cast<size_t>(w2.l.l*(w2.l.l + 1) + w2.l.m);

    if(at_i == at_j && l_i == l_j){
        res += XH1[at_i][l_i];
    }

    if(std::abs(w1.kappa*w1.kappa - w2.kappa*w2.kappa) < 1E-15){
        B1 = Bloch_summed_structure_constant(l_low_1, w1.kappa, cryst, w1.l, w2.l);
        res += -w2.kappa*w2.kappa*B1.dot_evaluate(tau_ij, kp);
    }

    B1 = Bloch_summed_structure_constant(l_low_2, w2.kappa, cryst, w1.l, w2.l);
    res += XH2[at_i][l_i]*B1.evaluate(tau_ij, kp);

    B1 = Bloch_summed_structure_constant(l_low_1, w1.kappa, cryst, w1.l, w2.l);
    res += B1.evaluate(tau_ij, kp)*XH2[at_j][l_j];
    if(std::abs(w1.kappa*w1.kappa - w2.kappa*w2.kappa) < 1E15){
        res += B1.evaluate(tau_ij, kp);
    }

    for(size_t at_m = 0; at_m < cryst.atoms.size(); at_m++){
        if(cryst.atoms[at_m].get_Z() > 20){
        }
        l_low = 3;
        for(Augmented_Bessel J1 : w1.J[at_m]){
            for(Augmented_Bessel J2 : w2.J[at_m]){
                if(J1.l == J2.l){
                    size_t ipp = static_cast<size_t>(J1.l.l*(J1.l.l + 1) + J1.l.m);
                    B1 = Bloch_summed_structure_constant(l_low, w1.kappa, cryst, J1.l, w1.l);
                    B2 = Bloch_summed_structure_constant(l_low, w2.kappa, cryst, J2.l, w2.l);
                    res += B1.evaluate(J1.center - w1.center.pos, kp).conjugate()*
                    XH3[at_m][ipp]*
                    B2.evaluate(J2.center - w2.center.pos, kp);
                }
            }
        }
    }
    return res;
}

double X_S1(const Augmented_Hankel& Ht1, const Augmented_Hankel& Ht2, const Atom& at)
{
    double res = 0;
    const Envelope_Hankel H1(at, Ht1.l, Ht1.kappa);
    const Envelope_Hankel H2(at, Ht2.l, Ht2.kappa);
    res += augmented_integral(Ht1, Ht2);
    res += off_atomic_integral(H1, H2);
    return res;
}

double X_S2(const Augmented_Hankel& Ht1, const Augmented_Bessel& Jt2, const Atom& at)
{
    double res = 0;
    const Envelope_Hankel H1(at, Ht1.l, Ht1.kappa);
    const Envelope_Bessel J2(at, Jt2.l, Jt2.kappa);
    // std::cout << "L = " << Ht1.l.l << " " << Ht1.l.m << std::endl;
    res += augmented_integral(Ht1, Jt2);
    // std::cout << "Augmented integral <Ht|Jt> = " << augmented_integral(Ht1, Jt2) << std::endl;
    res -= atomic_integral(H1, J2);
    // std::cout << "Atomic integral <H|J> = " << atomic_integral(H1, J2) << std::endl;
    if(Ht1.kappa*Ht1.kappa != Jt2.kappa*Jt2.kappa){
        res += 1./(-Jt2.kappa*Jt2.kappa + Ht1.kappa*Ht1.kappa);
    }
    // std::cout << "EJ = " << Jt2.EJ << " EH = " << Ht1.EH << std::endl;
    return res;
}

double X_S3(const Augmented_Bessel& Jt1, const Augmented_Bessel& Jt2, const Atom& at)
{
    double res = 0;
    const Envelope_Bessel J1(at, Jt1.l, Jt1.kappa);
    const Envelope_Bessel J2(at, Jt2.l, Jt2.kappa);
    // std::cout << "L = " << Jt1.l.l << " " << Jt1.l.m << std::endl;
    res += augmented_integral(Jt1, Jt2);
    // std::cout << "Augmented integral <Jt|Jt> = " << augmented_integral(Jt1, Jt2) << std::endl;
    res -= atomic_integral(J1, J2);
    // std::cout << "Atomic integral <J|J> = " << atomic_integral(J1, J2) << std::endl;
    return res;
}

GSL::Complex Simulation::S_element(const size_t i1, const size_t i2, const GSL::Vector& kp)
{
    const Augmented_spherical_wave w1 = basis_valence[i1];
    const Augmented_spherical_wave w2 = basis_valence[i2];
    Bloch_summed_structure_constant B1, B2;
    int l_low_1 = 2, l_low_2 = 2, l_low  = 2;
    if(w1.center.Z > 20){
        l_low_1 = 3;
    }
    if(w2.center.Z > 20){
        l_low_2 = 3;
    }
    GSL::Complex res = 0.;
    const GSL::Vector tau_ij(w1.center.pos - w2.center.pos);
    size_t at_i = 0, at_j = 0;
    for(size_t i = 0; i < cryst.atoms.size(); i++){
        if(cryst.atoms[i] == w1.center){
            at_i = i;
        }
        if(cryst.atoms[i] == w2.center){
            at_j = i;
        }
    }
    size_t l_i = static_cast<size_t>(w1.l.l*(w1.l.l + 1) + w1.l.m);
    size_t l_j = static_cast<size_t>(w2.l.l*(w2.l.l + 1) + w2.l.m);

    if(at_i == at_j && l_i == l_j){
        res += XS1[at_i][l_i];
    }

    if(std::abs(w1.kappa*w1.kappa - w2.kappa*w2.kappa) < 1E-15){
        B1 = Bloch_summed_structure_constant(l_low_1, w1.kappa, cryst, w1.l, w2.l);
        res += B1.dot_evaluate(tau_ij, kp);
    }

    B1 = Bloch_summed_structure_constant(l_low_1, w2.kappa, cryst, w1.l, w2.l);
    res += XS2[at_i][l_i]*B1.evaluate(tau_ij, kp);

    B1 = Bloch_summed_structure_constant(l_low_2, w1.kappa, cryst, w1.l, w2.l);
    res += B1.evaluate(tau_ij, kp)*XS2[at_j][l_j];

    for(size_t at_m = 0; at_m < cryst.atoms.size(); at_m++){
        if(cryst.atoms[at_m].get_Z() > 20){
            l_low = 3;
        }
        for(Augmented_Bessel J1 : w1.J[at_m]){
            for(Augmented_Bessel J2 : w2.J[at_m]){
                if(J1.l == J2.l){
                    size_t ipp = static_cast<size_t>(J1.l.l*(J1.l.l + 1) + J1.l.m);
                    B1 = Bloch_summed_structure_constant(l_low, w1.kappa, cryst, J1.l, w1.l);
                    B2 = Bloch_summed_structure_constant(l_low, w2.kappa, cryst, J2.l, w2.l);
                    res += B1.evaluate(J1.center - w1.center.pos, kp).conjugate()*
                    XS3[at_m][ipp]*
                    B2.evaluate(J2.center - w2.center.pos, kp);
                }
            }
        }
    }
    return res;
}


void Simulation::set_up_X_matrices()
{
    Augmented_spherical_wave w1, w2;
    size_t at_i = 0;

    std::cout << "Setting up X matrices (for H and S)" << std::endl;
    for(size_t i = 0; i < basis_valence.size(); i++){
        for(size_t j = 0; j < basis_valence.size(); j++){

            w1 = basis_valence[i];
            w2 = basis_valence[j];
            size_t l_i = static_cast<size_t>(w1.l.l*(w1.l.l + 1) + w1.l.m);

            for(size_t n_s = 0; n_s < cryst.atoms.size(); n_s++){
                if(cryst.atoms[n_s] == w1.center){
                    at_i = n_s;
                }
            }

            if(w1.center == w2.center && w1.H.l == w2.H.l){
                XH1[at_i][l_i] = X_H1(w1.H, w2.H, w1.center);
                XS1[at_i][l_i] = X_S1(w1.H, w2.H, w1.center);
            }
            // Find off center expansion on site w1.center of wave w2

            for(Augmented_Bessel J2s : w2.J[at_i]){
                if(J2s.l == w1.H.l){
                    XH2[at_i][l_i] = X_H2(w1.H, J2s, w1.center);
                    XS2[at_i][l_i] = X_S2(w1.H, J2s, w1.center);
                }
            }
            // Find off center expansion on site s of waves w1 ans w2
            for(size_t n_s = 0; n_s < cryst.atoms.size(); n_s++){
                for(Augmented_Bessel J1s : w1.J[n_s]){
                    for(Augmented_Bessel J2s : w2.J[n_s]){
                        if(J1s.l == J2s.l){
                            size_t ipp = static_cast<size_t>(J1s.l.l*(J1s.l.l + 1) + J1s.l.m);
                            XH3[n_s][ipp]
                              =  X_H3(J1s, J2s, cryst.atoms[n_s]);
                            XS3[n_s][ipp]
                              =  X_S3(J1s, J2s, cryst.atoms[n_s]);
                        }
                    }
                }
            }
        }
    }
    std::cout << "XS1 = " << XS1 << std::endl;
    std::cout << "XS2 = " << XS2 << std::endl;
    std::cout << "XS3 = " << XS3 << std::endl;
}

void Simulation::set_up_H(const GSL::Vector& kp)
{
    // std::cout << "Setting up Hamiltonian matrix." << std::endl;
    for(size_t i = 0; i < basis_valence.size(); i++){
        for(size_t j = 0; j < basis_valence.size(); j++){
            this->H[i].set(j, H_element(i, j, kp));
//            this->H[j].set(i, this->H[i][j].conjugate());
        }
    }
}

void Simulation::set_up_S(const GSL::Vector& kp)
{
    // std::cout << "Setting up overlap matrix." << std::endl;
    for(size_t i = 0; i < basis_valence.size(); i++){
        for(size_t j = 0; j < basis_valence.size(); j++){
            this->S[i].set(j, S_element(i, j, kp));
//            this->S[j].set(i, this->S[i][j].conjugate());
        }
    }
//    std::cout << "Overlap matrix" << std::endl;
//    for(size_t i = 0 ; i < basis_valence.size(); i++){
//        std::cout << "  " << S[i] << std::endl;
//    }
}

void Simulation::calc_eigen()
{
    size_t N = basis_valence.size()/2;
    GSL::Complex_Matrix eigvecs_up(N, N);
    GSL::Complex_Matrix eigvecs_down(N, N);
    GSL::Vector eigvals_up(N);
    GSL::Vector eigvals_down(N);

    GSL::Complex_Matrix H_up(N, N);
    GSL::Complex_Matrix H_down(N, N);
    GSL::Complex_Matrix S_up(N, N);
    GSL::Complex_Matrix S_down(N, N);

    for(size_t i = 0; i < N; i++){
        for(size_t j = 0; j < N; j++){
            H_up[i].set(j, H[2*i][2*j]);
            S_up[i].set(j, S[2*i][2*j]);
            H_down[i].set(j, H[2*i+1][2*j+1]);
            S_down[i].set(j, S[2*i+1][2*j+1]);
        }
    }

    std::pair<GSL::Complex_Matrix, GSL::Vector> tmp;
    tmp = GSL::hermitian_eigen(S_up);
//    std::cout << "Overlap matrix" << std::endl;
//    for(size_t i = 0 ; i < N; i++){
//        std::cout << "  " << S_up[i] << std::endl;
//    }
    std::cout << "Eigenvalues of S " << tmp.second << std::endl;
    tmp = GSL::general_hermitian_eigen(H_up, S_up);
    eigvecs_up = tmp.first;
    eigvals_up = tmp.second;
    tmp  = GSL::general_hermitian_eigen(H_down, S_down);
    eigvecs_down = tmp.first;
    eigvals_down = tmp.second;

    for(size_t i = 0; i < N; i++){
        std::cout << "Eigenvalues (up, down): " << eigvals_up[i] << " " <<
        eigvals_down[i] << "Ry\n";
        std::cout << "Eigenvectors (up, down) :\n";
        std::cout << "  " << eigvecs_up.get_col(i) << "\n";
        std::cout << "  " << eigvecs_down.get_col(i) << "\n";
    }
}
