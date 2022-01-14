#ifndef MIXER_H
#define MIXER_H

#include <vector>
#include <GSLpp/vector.h>

template<class AT_QUANTITY>
class Mixer {
protected:
    Mixer() = default;
    Mixer(const Mixer&) = default;
    Mixer(Mixer&&) = default;
public:
    virtual AT_QUANTITY mix(const AT_QUANTITY& in, const AT_QUANTITY& out) = 0;
};

template<class AT_QUANTITY>
class Linear_mixer : public Mixer<AT_QUANTITY> {
protected:
    double alpha_m;
    Linear_mixer() = default;
    Linear_mixer(const Linear_mixer&) = default;
    Linear_mixer(Linear_mixer&&) = default;
public:
    Linear_mixer(double a)
     : Mixer<AT_QUANTITY>(), alpha_m(a)
    {}

    virtual AT_QUANTITY mix(const AT_QUANTITY& in, const AT_QUANTITY& out);
};

template<>
class Linear_mixer<Density>: public Mixer<Density> {
protected:
    double alpha_m;
    Linear_mixer() = default;
    Linear_mixer(const Linear_mixer&) = default;
    Linear_mixer(Linear_mixer&&) = default;
public:
    Linear_mixer(double a)
     : Mixer<Density>(), alpha_m(a)
    {}

    virtual Density mix(const Density& in, const Density& out)
    {
        Density res({1001}, {0});
        for(size_t i = 0; i < res.n_at(); i++ ){
            res.valence(i) = out.valence(i)*alpha_m + in.valence(i)*(1 - alpha_m);
        }
        return res;
    }
};

template<class AT_QUANTITY>
class Linear_magnetic_mixer : public Mixer<AT_QUANTITY> {
protected:
    double alpha_m, beta_m;
    Linear_magnetic_mixer() = default;
    Linear_magnetic_mixer(const Linear_magnetic_mixer&) = default;
    Linear_magnetic_mixer(Linear_magnetic_mixer&&) = default;
public:
    Linear_magnetic_mixer(double a, double b)
     : Mixer<AT_QUANTITY>(), alpha_m(a), beta_m(a)
    {}

    virtual AT_QUANTITY mix(const AT_QUANTITY& in, const AT_QUANTITY& out);
};

template<>
class Linear_magnetic_mixer<Density>: public Mixer<Density> {
protected:
    double alpha_m, beta_m;
    Linear_magnetic_mixer() = default;
    Linear_magnetic_mixer(const Linear_magnetic_mixer&) = default;
    Linear_magnetic_mixer(Linear_magnetic_mixer&&) = default;
public:
    Linear_magnetic_mixer( double a, double b)
     : Mixer<Density>(), alpha_m(a), beta_m(b)
    {}

    virtual Density mix(const Density& in, const Density& out)
    {
        if(!in.spinpol()){
               throw std::runtime_error("Density must be spin polarized for magnetic mixer!");
        }
        Density res(in.mesh_lengths(), in.lmax(), true);
        for(size_t i = 0; i < res.n_at(); i++ ){
            size_t len = res.mesh_lengths()[i];
            for(size_t l = 0; l < res.lmax()[i] + 1; l++){
                    GSL::Vector::Const_View in_up(in.valence(i, l), 0, len, 2), in_down(in.valence(i, l), 1, len, 2),
                            out_up(out.valence(i, l), 0, len, 2), out_down(out.valence(i, l), 1, len, 2);
                    GSL::Vector::View res_up(res.valence(i, l), 0, len, 2), res_down(res.valence(i, l), 0, len, 2);

                    res_up = out_up*(alpha_m + beta_m) + out_down*(alpha_m - beta_m) + in_up*(2 - alpha_m - beta_m) + in_down*(beta_m - alpha_m);
                    res_down = out_up*(alpha_m - beta_m) + out_down*(alpha_m + beta_m) + in_up*(beta_m - alpha_m ) + in_down*(2 - alpha_m - beta_m);
                    res_up *= 0.5;
                    res_down *= 0.5;
            }
        }
        return res;
    }
};

#endif // MIXER_H
