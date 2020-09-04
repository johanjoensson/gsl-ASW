#ifndef EWALD_INT_H
#define EWALD_INT_H

#include <vector>
#include "log_mesh.h"
#include "utils.h"

class Ewald_integral{
        double bar_ew_int(const double kappa, const double eta, const lm l, const double r) const;
        double bar_comp_ew_int(const double kappa, const double eta, const lm l, const double r) const;
    public:
        Ewald_integral(){}
        std::vector<double> evaluate(const double kappa, const double eta, const lm l, const Logarithmic_mesh &mesh) const;
        std::vector<double> evaluate_comp(const double kappa, const double eta, const lm l, const Logarithmic_mesh &mesh) const;
        double ewald_int(const double kappa, const double eta, const lm l, const double r) const;
        double comp_ewald_int(const double kappa, const double eta, const lm l, const double r) const;

        double operator()(const double kappa, const double eta, const lm l, const double r) const
        {
            return ewald_int(kappa, eta, l, r);
        }
        double comp(const double kappa, const double eta, const lm l, const double r) const
        {
            return comp_ewald_int(kappa, eta, l, r);
        }
};

#endif //EWALD_INT_H
