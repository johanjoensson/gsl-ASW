#ifndef EWALD_INT_H
#define EWALD_INT_H

#include <vector>
#include "log_mesh.h"
#include "spherical_fun.h"

/***************************************************************************//**
* A class for the integral representation of Hankel functions\n
* Contains:\n
* __ewald_param__ - Parameter separating long range part from short range one\n
* __kappa__ - Energy parameter used (usually kappa^2 = -0.015)\n
* __h__ -
*******************************************************************************/
class Ewald_integral{
        double ewald_param;
        double kappa;
        double bar_ew_int(lm l, double r);
        double bar_comp_ew_int(lm l, double r);
    public:
        //! Set the Ewald parameter
        void set_ewald_param(double eta);
        //! Set energy parameter
        void set_kappa(double kappa);

        Ewald_integral();

        double ewald_int(lm l, double r);
        double comp_ewald_int(lm l, double r);
        //! Calculate the Ewald integral for angular quantum number (l.l, l.m) and evaluate it on every point in the mesh
        std::vector<double> evaluate(lm l, Logarithmic_mesh &mesh);
        //! Calculate the complementary Ewald integral for angular quantum number (l.l, l.m) and evaluate it on every point in the mesh
        std::vector<double> evaluate_comp(lm l, Logarithmic_mesh &mesh);
        double ewald_int(lm l, double r);
        double comp_ewald_int(lm l, double r);
};

#endif //EWALD_INT_H
