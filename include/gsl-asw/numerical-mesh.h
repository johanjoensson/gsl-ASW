#ifndef NUMERICAL_MESH_LIB_H
#define NUMERICAL_MESH_LIB_H
#include <array>
#include <vector>
#include <cmath>
#include <algorithm>
#include <type_traits>
#include <tuple>
#include <iostream>

#if __cplusplus >= 202002L
concept arithmetic = std::integral || std::floating_point
#endif

namespace{
#if __cplusplus >= 202002L
    template<arithmetic T>
#else
    template<typename T>
#endif
    constexpr int sgn(T val)
    {
    return (T(0) < val) - (val < T(0));
    }
}

#if __cplusplus >= 202002L
template<size_t Dim, floating_point Scalar = double>
#else
template<size_t Dim, class Scalar = double>
#endif
class Mesh_base {
    using Arr = std::array<Scalar, Dim>;

protected:
    template<class A, class B>
    std::array<std::tuple<A, B>, Dim> zip2 (const std::array<A, Dim>& a, const std::array<B, Dim>& b) const
    {
        std::array<std::tuple<A, B>, Dim> res;
#if __cplusplus >= 202002L
        std::transform(std::begin(a), std::end(a), std::begin(b), std::begin(res),
            [](const A& ai, const B& bi)
            {
                return std::tuple<A, B>(ai, bi);
            });
#else
        auto ai = a.begin();
        auto bi = b.begin();
        auto ri = res.begin();
        for(;ai != a.end(); ai++, bi++, ri++){
            *ri = {*ai, *bi};
        }
#endif
            return res;
    }

    template<class A, class B, class C>
    std::array<std::tuple<A, B, C>, Dim> zip3(const std::array<A, Dim>& a, const std::array<B, Dim>& b, std::array<C, Dim> c) const
    {
        auto tmp = zip2(b, c);
        std::array<std::tuple<A, B, C>, Dim> res;
#if __cplusplus >= 202002L
        std::transform(std::begin(a), std::end(a), std::begin(tmp), std::begin(res),
            [] (const Scalar ai, const std::tuple<B, C> ti)
            {
                auto [bi, ci] = ti;
                return std::tuple<A, B, C>(ai, bi, ci);
            });
#else
        auto ai = a.begin();
        auto bi = b.begin();
        auto ci = c.begin();
        auto ri = res.begin();
        for(;ai != a.end(); ai++, bi++, ci++, ri++){
            *ri = {*ai, *bi, *ci};
        }
#endif
        return res;
    }

    template<class A, class B, class C, class D>
    std::array<std::tuple<A, B, C, D>, Dim> zip4 (const std::array<A, Dim>& a, const std::array<B, Dim>& b, const std::array<C, Dim>& c, const std::array<D, Dim>& d) const
    {
        auto z1 = zip2(a, b), z2 = zip2(c, d);
        std::array<std::tuple<A, B, C, D>, Dim> res;
#if __cplusplus >= 202002L
        std::transform(std::begin(z1), std::end(z1), std::begin(z2), std::begin(res),
        [] (const std::tuple<A, B> ti, const std::tuple<C, D> ui)
        {
            auto [ai, bi] = ti;
            auto [ci, di] = ui;
            return std::tuple<A, B, C, D>(ai, bi, ci, di);
        });
#else
        auto ai = a.begin();
        auto bi = b.begin();
        auto ci = c.begin();
        auto di = d.begin();
        auto ri = res.begin();
        for(;ai != a.end(); ai++, bi++, ci++, di++, ri++){
            *ri = {*ai, *bi, *ci, *di};
        }
#endif
        return res;
    }

    template<class A, class B, class C, class D>
    std::array<std::tuple<A, B, C, D>, Dim> zip4 (const std::array<A, Dim>& a, const std::array<std::tuple<B, C, D>, Dim>& t) const
    {
        std::array<std::tuple<A, B, C, D>, Dim> res;
#if __cplusplus >= 202002L
        std::transform(std::begin(a), std::end(a), std::begin(t), std::begin(res),
        [] (const A ai, const std::tuple<B, C, D> ti)
        {
            auto [bi, ci, di] = ti;
            return std::tuple<A, B, C, D>(ai, bi, ci, di);
        });
#else
        auto ai = a.begin();
        auto ti = t.begin();
        auto ri = res.begin();
        for(;ai != a.end(); ai++, ti++, ri++){
            auto bi = std::get<0>(*ti);
            auto ci = std::get<1>(*ti);
            auto di = std::get<2>(*ti);
            *ri = {*ai, bi, ci, di};
        }
#endif
        return res;
    }

    Arr alpha_m, beta_m, gamma_m;
    std::array<size_t, Dim> N_m;
    Mesh_base() = default;

public:
    virtual ~Mesh_base() = default;
    Mesh_base(const Arr alpha, const Arr beta, const Arr gamma, const std::array<size_t, Dim>& N)
     : alpha_m(alpha), beta_m(beta), gamma_m(gamma), N_m(N)
    {}

    Mesh_base(const Mesh_base&) = default;
    Mesh_base(Mesh_base&&) = default;

    Mesh_base& operator=(const Mesh_base&) = default;
    Mesh_base& operator=(Mesh_base&&) = default;

    auto parameters() const {return zip3(alpha_m, beta_m, gamma_m);}

    virtual Arr r(const Arr& x) const = 0;
    virtual Arr r2(const Arr& x) const = 0;
    virtual Arr dr(const Arr& x) const = 0;

    Arr operator()(const Arr& x) const {return r(x);}


    std::array<std::vector<Scalar>, Dim> r() const;

    std::array<size_t, Dim> dim() const {return N_m;}

};

template<class Scalar>
class Mesh_base<1, Scalar> {
protected:
    Scalar alpha_m, beta_m, gamma_m;
    size_t N_m;
    Mesh_base() = default;

public:
    virtual ~Mesh_base() = default;
    Mesh_base(const Scalar alpha, const Scalar beta, const Scalar gamma, const size_t N)
     : alpha_m(alpha), beta_m(beta), gamma_m(gamma), N_m(N)
    {}

    Mesh_base(const Mesh_base&) = default;
    Mesh_base(Mesh_base&&) = default;

    Mesh_base& operator=(const Mesh_base&) = default;
    Mesh_base& operator=(Mesh_base&&) = default;

    virtual Scalar r(const Scalar& x) const = 0;
    virtual Scalar r2(const Scalar& x) const = 0;
    virtual Scalar dr(const Scalar& x) const = 0;
    Scalar operator()(const Scalar& x) const {return r(x);}

    virtual std::vector<Scalar> r() const
    {
        std::vector<Scalar> res;
        for(size_t i = 0; i < this->N_m; i++){
            res.push_back(this->r(i));
        }
        return res;
    }

    virtual std::vector<Scalar> r2() const
    {
        std::vector<Scalar> res;
        for(size_t i = 0; i < this->N_m; i++){
            res.push_back(this->r2(i));
        }
        return res;
    }

    virtual std::vector<Scalar> dr() const
    {
        std::vector<Scalar> res;
        for(size_t i = 0; i < this->N_m; i++){
            res.push_back(this->dr(i));
        }
        return res;
    }

    template<class Func>
    std::vector<Scalar> evaluate(Func f, size_t h = 1)
    {
        std::vector<Scalar> res;
        for(size_t i = 0; i < this->N_m; i += h){
            res.push_back(f(this->r(i)));
        }
        return res;
    }

    template<class Func>
    std::vector<Scalar> evaluate(Func f, Scalar x0, Scalar x1, Scalar h)
    {
        std::vector<Scalar> res;
        for(Scalar x = x0; x <= x1; x += h){
            res.push_back(f(this->r(x)));
        }
        return res;
    }

    size_t dim() const {return N_m;}
};

template<size_t D, class Scalar = double>
class Linear_mesh : public Mesh_base<D, Scalar> {
private:
    Linear_mesh() = default;

public:
    using Arr = std::array<Scalar, D>;
    Linear_mesh(const Arr& R_min, const Arr& R_max, const std::array<size_t, D>& N)
    : Mesh_base<D, Scalar>(
        [this] (const Arr& R_min, const Arr& R_max,  const std::array<size_t, D>& N)
        {
            Arr res;
            auto z3 = this->zip3(R_min,  R_max, N);
            std::transform(std::begin(z3), std::end(z3), std::begin(res),
                [] (std::tuple<Scalar, Scalar, size_t> ti)
                {
#if __cplusplus >= 202002L
                    auto [r1, r2, n] = ti;
#else
                    auto r1 = std::get<0>(ti);
                    auto r2 = std::get<1>(ti);
                    auto n = std::get<2>(ti);
#endif
                    return (r2 - r1)/(n - 1);
                });
                return res;
            }(R_min, R_max, N)
        , Arr(), R_min, N)
    {}

    ~Linear_mesh() = default;

    Linear_mesh(const Linear_mesh&) = default;
    Linear_mesh(Linear_mesh&&) = default;

    Linear_mesh& operator=(const Linear_mesh&) = default;
    Linear_mesh& operator=(Linear_mesh&&) = default;

    Arr r(const Arr& x) const override
    {
        Arr res;

#if __cplusplus >= 202002L
        std::transform(std::begin(x), std::end(x), std::begin(this->parameters()), std::begin(res),
        [] (const Scalar xi, const std::tuple<Scalar, Scalar, Scalar> ti)
        {
            auto [ai, bi, ci] = ti;
            return ai*xi +0*bi + ci;
        });
#else
        auto z4 = this->zip4(x, this->parameters());
        std::transform(std::begin(z4), std::end(z4), std::begin(res),
        [] (const std::tuple<Scalar, Scalar, Scalar, Scalar> ti)
        {
            auto xi = std::get<0>(ti);
            auto ai = std::get<1>(ti);
            auto bi = std::get<2>(ti);
            auto ci = std::get<3>(ti);
            return ai*xi +0*bi + ci;
        });
#endif
        return res;
    }

    Arr r2(const Arr& x) const override
    {
        Arr res, r = this->r(x);
        std::transform(std::begin(r), std::end(r), std::begin(res),
            [] (const Scalar s) {return s*s;});
        return res;
    }

    Arr dr(const Arr& x) const override
    {
        return this->alpha_m;
    }
};

template<class Scalar>
class Linear_mesh<1, Scalar> : public Mesh_base<1, Scalar> {
private:
    Linear_mesh() = default;

public:
    Linear_mesh(const Scalar R_min, const Scalar R_max, const size_t N)
    : Mesh_base<1, Scalar>((R_max - R_min)/(N - 1), 0, R_min, N)
    {}

    ~Linear_mesh() = default;

    Linear_mesh(const Linear_mesh&) = default;
    Linear_mesh(Linear_mesh&&) = default;

    Linear_mesh& operator=(const Linear_mesh&) = default;
    Linear_mesh& operator=(Linear_mesh&&) = default;

    using Mesh_base<1, Scalar>::r;
    using Mesh_base<1, Scalar>::r2;
    using Mesh_base<1, Scalar>::dr;

    Scalar r(const Scalar& x) const override
    {
        return x*this->alpha_m + this->gamma_m;
    }

    Scalar r2(const Scalar& x) const override
    {
        return std::pow(this->r(x), 2);
    }

    Scalar dr(const Scalar& x) const override
    {
        return this->alpha_m;
    }
};

template<size_t D, class Scalar = double>
class Quadratic_mesh : public Mesh_base<D, Scalar> {
public:
    using Arr = std::array<Scalar, D>;
    Quadratic_mesh() = delete;
    Quadratic_mesh(const Arr& R_min, const Arr& R_max, const std::array<size_t, D>& N)
    : Mesh_base<D, Scalar>(
        [this] (const Arr& R_min, const Arr& R_max, const std::array<size_t, D>& N)
        {
            Arr res;
            auto z3 = this->zip3(R_min, R_max, N);
            std::transform(std::begin(z3), std::end(z3), std::begin(res),
                [] (const std::tuple<Scalar, Scalar, size_t>& ti)
                {
#if __cplusplus >= 202002L
                    auto [r1, r2, n] = ti;
#else
                    auto r1 = std::get<0>(ti);
                    auto r2 = std::get<1>(ti);
                    auto n = std::get<2>(ti);
#endif
                    return (r2-r1)/std::pow(n - 1, 2);
                });
            return res;
        }(R_min, R_max, N)
        , Arr(), R_min, N )
    {}

    ~Quadratic_mesh() = default;

    Quadratic_mesh(const Quadratic_mesh&) = default;
    Quadratic_mesh(Quadratic_mesh&&) = default;

    Quadratic_mesh& operator=(const Quadratic_mesh&) = default;
    Quadratic_mesh& operator=(Quadratic_mesh&&) = default;

    Arr r(const Arr& x) const override
    {
        Arr res;
#if __cplusplus >= 202002L
        std::transform(std::begin(x), std::end(x), std::begin(this->parameters()), std::begin(res),
            [] (const Scalar& xi, const std::tuple<Scalar, Scalar, Scalar>& ti)
            {
                auto [a, b, c] = ti;
                return sgn(xi)*a*std::pow(xi, 2) + c;
            });
#else
        auto z4 = this->zip4(x, this->parameters());
        std::transform(std::begin(z4), std::end(z4), std::begin(res),
            [] (const std::tuple<Scalar, Scalar, Scalar, Scalar>& ti)
            {
                auto xi = std::get<0>(ti);
                auto a = std::get<1>(ti);
                auto b = std::get<2>(ti);
                auto c = std::get<3>(ti);
                return sgn(xi)*a*std::pow(xi, 2) + c;
            });
#endif
        return res;
    }

    Arr r2(const Arr& x) const override
    {
        Arr res, r = this->r(x);
        std::transform(std::begin(r), std::end(r), std::begin(res), [] (const Scalar& s) {return s*s;});
        return res;
    }

    Arr dr(const Arr& x) const override
    {
        Arr res;
#if __cplusplus >= 202002L
        std::transform(std::begin(x), std::end(x), std::begin(this->parameters()), std::begin(res),
            [] (const Scalar xi, const std::tuple<Scalar, Scalar, Scalar> ti)
            {
                auto [a, b, c] = ti;
                return sgn(xi)*2*a*xi + b*0 + c*0;
            });
#else
        auto z4 = this->zip4(x, this->parameters());
        std::transform(std::begin(z4), std::end(z4), std::begin(res),
            [] (const std::tuple<Scalar, Scalar, Scalar, Scalar>& ti)
            {
                auto xi = std::get<0>(ti);
                auto a = std::get<1>(ti);
                auto b = std::get<2>(ti);
                auto c = std::get<3>(ti);
                return sgn(xi)*2*a*xi + b*0 + c*0;
            });
#endif
        return res;
    }
};

template<class Scalar>
class Quadratic_mesh<1, Scalar> : public Mesh_base<1, Scalar> {
public:
    Quadratic_mesh() = delete;
    Quadratic_mesh(const Scalar R_min, const Scalar R_max, const size_t N)
    : Mesh_base<1, Scalar>((R_max - R_min)/std::pow(N - 1, 2), 0, R_min, N )
    {}

    ~Quadratic_mesh() = default;

    Quadratic_mesh(const Quadratic_mesh&) = default;
    Quadratic_mesh(Quadratic_mesh&&) = default;

    Quadratic_mesh& operator=(const Quadratic_mesh&) = default;
    Quadratic_mesh& operator=(Quadratic_mesh&&) = default;

    Scalar r(const Scalar& x) const override
    {
        return sgn(x)*std::pow(x, 2)*this->alpha_m + this->gamma_m;
    }

    Scalar r2(const Scalar& x) const override
    {
        return std::pow(this->r(x), 2);
    }

    Scalar dr(const Scalar& x) const override
    {
        return 2*this->alpha_m*sgn(x)*x;
    }
};

template<size_t D, class Scalar = double>
class Exponential_mesh : public Mesh_base<D, Scalar> {
public:
    using Arr = std::array<Scalar, D>;
    Exponential_mesh() = delete;
    Exponential_mesh(const Arr& R_min, const Arr& R_max, const Arr& beta, const std::array<size_t, D>& N)
    : Mesh_base<D, Scalar>([this] (const Arr& R_min, const Arr& R_max, const Arr& beta, const std::array<size_t, D>& N)
        {
            Arr res;
            auto z4 = this->zip4(R_min, R_max, beta, N);
            std::transform(std::begin(z4), std::end(z4), std::begin(res),
                [] (const std::tuple<Scalar, Scalar, Scalar, size_t> ti)
                {
#if __cplusplus >= 202002L
                    auto [r1, r2, b, n] = ti;
#else
                    auto r1 = std::get<0>(ti);
                    auto r2 = std::get<1>(ti);
                    auto b = std::get<2>(ti);
                    auto n = std::get<2>(ti);
#endif
                    return (r2 - r1)/(std::exp(b*(n - 1)) - 1);
                });
            return res;
        }(R_min, R_max, beta, N), beta, R_min, N)
    {}

    ~Exponential_mesh() = default;

    Exponential_mesh(const Exponential_mesh&) = default;
    Exponential_mesh(Exponential_mesh&&) = default;

    Exponential_mesh& operator=(const Exponential_mesh&) = default;
    Exponential_mesh& operator=(Exponential_mesh&&) = default;

    Arr r(const Arr& x) const override
    {
        Arr res;
#if __cplusplus >= 202002L
        std::transform(std::begin(x), std::end(x), std::begin(this->parameters()), std::begin(res),
            [] (const Scalar& xi, const std::tuple<Scalar, Scalar, Scalar>& ti)
            {
                auto [a, b, c] = ti;
                return sgn(xi)*a*(std::exp(sgn(xi)*xi*b) - 1) + c;
            });
#else
        auto z4 = this->zip4(x, this->parameters());
        std::transform(std::begin(z4), std::end(z4), std::begin(res),
            [] (const std::tuple<Scalar, Scalar, Scalar, Scalar>& ti)
            {
                auto xi = std::get<0>(ti);
                auto a = std::get<1>(ti);
                auto b = std::get<2>(ti);
                auto c = std::get<3>(ti);
                return sgn(xi)*a*(std::exp(sgn(xi)*xi*b) - 1) + c;
            });
#endif
        return res;
    }

    Arr r2(const Arr& x) const override
    {
        Arr res, r = this->r(x);
        std::transform(std::begin(r), std::end(r), std::begin(res), [] (const Scalar& s) {return s*s;});
        return res;
    }

    Arr dr(const Arr& x) const override
    {
        Arr x_abs;
        std::transform(std::begin(x), std::end(x), std::begin(x_abs),
            [] (const Scalar& xi)
            {
                return sgn(xi)*xi;
            }
        );
        Arr r = this->r(x_abs);

        Arr res;
#if __cplusplus >= 202002L
        std::transform(std::begin(r), std::end(r), std::begin(this->parameters()), std::begin(res),
            [] (const Scalar& ri, const std::tuple<Scalar, Scalar, Scalar>& ti)
            {
                auto [a, b, c] = ti;
                return b*(ri + a - c);
            });
#else
        auto z4 = this->zip4(r, this->parameters());
        std::transform(std::begin(z4), std::end(z4), std::begin(res),
            [] (const std::tuple<Scalar, Scalar, Scalar, Scalar>& ti)
            {
                auto ri = std::get<0>(ti);
                auto a = std::get<1>(ti);
                auto b = std::get<2>(ti);
                auto c = std::get<3>(ti);
                return b*(ri + a - c);
            });
#endif
        return res;
    }
};

template<class Scalar>
class Exponential_mesh<1, Scalar> : public Mesh_base<1, Scalar> {
public:
    Exponential_mesh() = delete;
    Exponential_mesh(const Scalar R_min, const Scalar R_max, const Scalar beta, const size_t N)
    : Mesh_base<1, Scalar>((R_max - R_min)/std::expm1(beta*(N-1)), beta, R_min, N)
    {}

    ~Exponential_mesh() = default;

    Exponential_mesh(const Exponential_mesh&) = default;
    Exponential_mesh(Exponential_mesh&&) = default;

    Exponential_mesh& operator=(const Exponential_mesh&) = default;
    Exponential_mesh& operator=(Exponential_mesh&&) = default;

    Scalar r(const Scalar& x) const override
    {
        Scalar a = this->alpha_m, b = this->beta_m, c = this->gamma_m;
        return sgn(x)*a*std::expm1(sgn(x)*b*x) + c;
    }

    Scalar r2(const Scalar& x) const override
    {
        return std::pow(this->r(x), 2);
    }

    Scalar dr(const Scalar& x) const override
    {
        Scalar a = this->alpha_m, b = this->beta_m, c = this->gamma_m;
        return b*(this->r(sgn(x)*x) + a - c);
    }
};

#endif //NUMERICAL_MESH_LIB_H
