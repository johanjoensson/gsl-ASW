#ifndef EWALD_INT_H
#define EWALD_INT_H

#include <vector>
#include "log_mesh.h"
#include "utils.h"

class Ewald_integral{
        double ewald_param;
        double kappa;
        double bar_ew_int(const lm l, const double r) const;
        double bar_comp_ew_int(const lm l, const double r) const;
    public:

        void set_ewald_param(const double eta);
        void set_kappa(const double kappa);

        Ewald_integral();
        std::vector<double> evaluate(const lm l, const Logarithmic_mesh &mesh) const;
        std::vector<double> evaluate_comp(const lm l, const Logarithmic_mesh &mesh) const;
        double ewald_int(const lm l, const double r) const;
        double comp_ewald_int(const lm l, const double r) const;
};

#endif //EWALD_INT_H
