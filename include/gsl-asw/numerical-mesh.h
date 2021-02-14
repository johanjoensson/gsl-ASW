#ifndef NUMERICAL_MESH_LIB_H
#define NUMERICAL_MESH_LIB_H

#if __cplusplus < 201103L // Require at least c++11
#error The numerical-mesh library requires at least c++11!
#endif

#include <array>
#include <vector>
#include <cmath>
#include <algorithm>
#include <type_traits>
#include <tuple>
#include <iostream>

/*******************************************************************************
* The __cplusplus macro is used to determine which c++ standard is used for the*
* compilation. This is "needed" because I want to use this library in projects *
* using c++11 as well as projects using c++20. This means parts of the code had*
* to be rewritten in order to enable the cool c++20 features only if c++20 is  *
* used. As a side effect this also made it possible to write code for c++17 as *
* well.                                                                        *
* Features requiring c++17:                                                    *
*  -- structured bindings i.e., auto [a,b,c,...] = ...                         *
* Features requiring c++20:                                                    *
*  -- Concepts e.g., integral and floating_point                               *
*  -- ranges-versions of algorithms                                            *
*  -- std::transform with 2 input containers                                   *
*******************************************************************************/

/*!*****************************************************************************
 * \defgroup NumMesh Numerical mesh
 * @{
 ******************************************************************************/

/*!*****************************************************************************
* The arithmetic concept, integral or floating point, is used to signify any
* numeric type.
*******************************************************************************/
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

/*!*****************************************************************************
* Base class for all numerical mesh classes. All meshes can be described by
* (max) three parameters and the number of mesh points (including endpoints)
* for each direction.
*******************************************************************************/
#if __cplusplus >= 202002L
template<size_t Dim, floating_point Scalar = double>
#else
template<size_t Dim, class Scalar = double>
#endif
class Mesh_base {
    using Arr = std::array<Scalar, Dim>;

    void increment(Arr& arr, Scalar step)
        const noexcept
    {
        arr[0] += step;
        for(auto j = 0; j < Dim - 1; j++){
            while(arr[j] >= N_m[j]){
                arr[j + 1]++;
                arr[j] -= N_m[j];
            }
        }
    }

    void decrement(Arr& arr, Scalar step)
        const noexcept
    {
        arr[0] -= step;
        for(auto j = 0; j < Dim; j++){
            while(arr[j] < 0){
                arr[j + 1]--;
                arr[j] += N_m[j];
            }
        }
    }

protected:
    template<class A, class B>
    std::array<std::tuple<A, B>, Dim> zip2 (const std::array<A, Dim>& a, const std::array<B, Dim>& b)
        const noexcept
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
    std::array<std::tuple<A, B, C>, Dim> zip3(const std::array<A, Dim>& a,
        const std::array<B, Dim>& b, std::array<C, Dim> c)
        const noexcept
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
    std::array<std::tuple<A, B, C, D>, Dim> zip4 (const std::array<A, Dim>& a,
        const std::array<B, Dim>& b, const std::array<C, Dim>& c, const std::array<D, Dim>& d)
        const noexcept
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
    std::array<std::tuple<A, B, C, D>, Dim> zip4 (const std::array<A, Dim>& a,
        const std::array<std::tuple<B, C, D>, Dim>& t)
        const noexcept
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

    //! Constructor
    /*!*************************************************************************
    * Construct a mesh from three arrays of parameters, and an array for the
    * number of mesh points along the different directions.
    * @param alpha First parameter for the mesh.
    * @param beta Second parameter for the mesh.
    * @param gamma Third parameter for the mesh (usually the end points of the
    * mesh).
    * @param N Number of mesh points in the mesh (including endpoints).
    ***************************************************************************/
    Mesh_base(const Arr alpha, const Arr beta, const Arr gamma, const std::array<size_t, Dim>& N)
     noexcept
     : alpha_m(alpha), beta_m(beta), gamma_m(gamma), N_m(N)
    {}

    Mesh_base(const Mesh_base&) = default;
    Mesh_base(Mesh_base&&) = default;

    Mesh_base& operator=(const Mesh_base&) = default;
    Mesh_base& operator=(Mesh_base&&) = default;

    std::array<std::tuple<Scalar, Scalar, Scalar>, Dim> parameters() const noexcept {return zip3(alpha_m, beta_m, gamma_m);}


public:
    virtual ~Mesh_base() = default;

    /*!*************************************************************************
    * Virtual function returning position of mesh point at x.
    * To be overwritten in subclasses.
    ***************************************************************************/
    virtual Arr r(const Arr& x) const noexcept = 0;

    /*!*************************************************************************
    * Virtual function returning squared position of mesh point at x.
    * To be overwritten in subclasses.
    ***************************************************************************/
    virtual Arr r2(const Arr& x) const noexcept = 0;

    /*!*************************************************************************
    * Virtual function returning step size of mesh at x.
    * To be overwritten in subclasses.
    ***************************************************************************/
    virtual Arr dr(const Arr& x) const noexcept = 0;

    Arr operator()(const Arr& x) const noexcept {return r(x);}

    //! Get all positions
    /*!*************************************************************************
    * Returns a vector containing all positions in the mesh.
    ***************************************************************************/
    std::vector<Arr> r()
        const noexcept
    {
        std::vector<Arr> res;
        Arr coord;
        coord.fill(0);
        while(coord.back() <  N_m.back()){
            res.push_back( r(coord) );
            increment(coord, 1);
        }
        return res;
    }

    //! Get mesh dimensions
    /*!*************************************************************************
    * Returns an array containing the number of mesh points along each directon.
    ***************************************************************************/
    std::array<size_t, Dim> dim() const noexcept { return N_m;}

    /*!*************************************************************************
    * Class for representing a generic point in a mesh.
    ***************************************************************************/
    class mesh_point {
    private:
        const Mesh_base* m_m;
        Arr i_m;
        mesh_point() = default;
    public:
        mesh_point(const Mesh_base* m, const Arr i) : m_m{m}, i_m{i}{}
        mesh_point(const mesh_point&) = default;
        mesh_point(mesh_point&&) = default;
        ~mesh_point() = default;

        mesh_point& operator=(const mesh_point&) = default;
        mesh_point& operator=(mesh_point&&) = default;

        Arr r() const noexcept {return m_m->r(i_m);}
        Arr r2() const noexcept {return m_m->r2(i_m);}
        Arr dr() const noexcept {return m_m->dr(i_m);}
    };


    /*!*************************************************************************
    * Class for enabling iteration through a mesh, constant version
    ***************************************************************************/
    class const_iterator{
    private:
        const Mesh_base* m_m;
        Arr i_m;
        const_iterator() = default;
    public:
        typedef const Arr difference_type;
        typedef const mesh_point value_type;
        // typedef const mesh_point reference;
        // typedef const mesh_point pointer;
        // typedef std::random_access_iterator_tag iterator_category;

        const_iterator(const const_iterator&) = default;
        const_iterator(const_iterator&&) = default;
        const_iterator(const Mesh_base* m, const Arr i) : m_m(m), i_m(i){}
        ~const_iterator() = default;
        const_iterator& operator=(const const_iterator&) = default;
        const_iterator& operator=(const_iterator&&) = default;

        value_type operator*() noexcept {return mesh_point(m_m, i_m);}

        bool operator==(const const_iterator& b)const noexcept{return i_m == b.i_m;}
        bool operator!=(const const_iterator& b)const noexcept{return !(*this == b);}
        bool operator>=(const const_iterator& b)const noexcept{return i_m >= b.i_m;}
        bool operator<=(const const_iterator& b)const noexcept{return i_m <= b.i_m;}
        bool operator>(const const_iterator& b)const noexcept{return !(*this <= b);}
        bool operator<(const const_iterator& b)const noexcept{return !(*this >= b);}


        const_iterator operator++(int)
            noexcept
        {
            auto res = *this;
            m_m->increment(i_m, 1);
            return res;
        }

        const_iterator& operator++()
            noexcept
        {
            m_m->increment(i_m, 1);
            return *this;
        }

        const_iterator operator--(int)
            noexcept
        {
            auto res = *this;
            m_m->decrement(i_m, 1);
            return res;
        }

        const_iterator& operator--()
            noexcept
        {
            m_m->decrement(i_m, 1);
            return *this;
        }

        const_iterator& operator+=(const difference_type n)
            noexcept
        {
            m_m->increment(i_m, n);
            return *this;
        }

        const_iterator& operator-=(const difference_type n)
            noexcept
        {
            m_m->decrement(i_m, n);
            return *this;
        }

        const_iterator operator+(const difference_type n)const noexcept{auto res = *this; return (res += n);}
        const_iterator operator-(const difference_type n)const noexcept{auto res = *this; return (res -= n);}
        difference_type operator-(const const_iterator& it)
            const noexcept
        {
            size_t diff = 0;
            Scalar offset = 1;
            for(auto j = 0; j < Dim; j++){
                diff += (i_m[j] - it.i_m[j])*offset;
                offset = 1;
                for(auto k = 0; k < j; k++){
                    offset *= N_m[k];
                }
            }
            return this->i_m - it.i_m;
        }

    };

    /*!*************************************************************************
    * Class for enabling reverse iteration through a mesh, constant version
    ***************************************************************************/
    typedef std::reverse_iterator<const_iterator> const_reverse_iterator;

    const_iterator begin()
        const noexcept
    {
        Arr i;
        i.fill(0);
        return const_iterator(this, i);
    }

    const_iterator cbegin()
        const noexcept
    {
        Arr i;
        i.fill(0);
        return const_iterator(this, i);
    }

    const_iterator end()
        const noexcept
    {
        Arr i;
        i.fill(0);
        i.back() = N_m.back();
        return const_iterator(this, i);
    }

    const_iterator cend()
        const  noexcept
    {
        Arr i;
        i.fill(0);
        i.back() = N_m.back();
        return const_iterator(this, i);
    }

    // reverse_iterator rbegin() noexcept{return reverse_iterator(this->end());}
    const_reverse_iterator rbegin() const noexcept{return const_reverse_iterator(this->end());}
    const_reverse_iterator crbegin() const noexcept{return const_reverse_iterator(this->cend());}
    // reverse_iterator rend() noexcept{return reverse_iterator(this->begin());}
    const_reverse_iterator rend() const noexcept{return const_reverse_iterator(this->begin());}
    const_reverse_iterator crend() const noexcept{return const_reverse_iterator(this->cbegin());}

};

/*!*****************************************************************************
* 1D template specialization, instead of using arrays we can use scalar values.
*******************************************************************************/
#if __cplusplus >= 202002L
template<floating_point Scalar>
#else
template<class Scalar>
#endif
class Mesh_base<1, Scalar> {
protected:
    Scalar alpha_m, beta_m, gamma_m;
    size_t N_m;
    Mesh_base() = default;

public:
    virtual ~Mesh_base() = default;
    //! Constructor
    /*!*************************************************************************
    * Construct a mesh from three parameters, and the number of mesh points.
    * @param alpha First parameter for the mesh.
    * @param beta Second parameter for the mesh.
    * @param gamma Third parameter for the mesh.
    * @param N Number of mesh points (including endpoints).
    ***************************************************************************/
    Mesh_base(const Scalar alpha, const Scalar beta, const Scalar gamma, const size_t N)
     noexcept
     : alpha_m(alpha), beta_m(beta), gamma_m(gamma), N_m(N)
    {}

    Mesh_base(const Mesh_base&) = default;
    Mesh_base(Mesh_base&&) = default;

    Mesh_base& operator=(const Mesh_base&) = default;
    Mesh_base& operator=(Mesh_base&&) = default;

    /*!*************************************************************************
    * Virtual function returning position of mesh point at x.
    * To be overwritten in subclasses.
    ***************************************************************************/
    virtual Scalar r(const Scalar& x) const noexcept = 0;

    /*!*************************************************************************
    * Virtual function returning square of position of mesh point at x.
    * To be overwritten in subclasses.
    ***************************************************************************/
    virtual Scalar r2(const Scalar& x) const noexcept = 0;

    /*!*************************************************************************
    * Virtual function returning step sze of mesh at x.
    * To be overwritten in subclasses.
    ***************************************************************************/
    virtual Scalar dr(const Scalar& x) const noexcept = 0;

    Scalar operator()(const Scalar& x) const noexcept {return r(x);}

    //! Get all positions
    /*!*************************************************************************
    * Returns a vector containing all positions in the mesh.
    ***************************************************************************/
    std::vector<Scalar> r()
        const noexcept
    {
        std::vector<Scalar> res;
        for(size_t i = 0; i < this->N_m; i++){
            res.push_back(this->r(i));
        }
        return res;
    }

    //! Get all positions squared
    /*!*************************************************************************
    * Returns a vector containing all squared positions in the mesh.
    ***************************************************************************/
    std::vector<Scalar> r2()
        const noexcept
    {
        std::vector<Scalar> res;
        for(size_t i = 0; i < this->N_m; i++){
            res.push_back(this->r2(i));
        }
        return res;
    }

    //! Get all step sizes
    /*!*************************************************************************
    * Returns a vector containing all step sizes in the mesh.
    ***************************************************************************/
    std::vector<Scalar> dr()
        const noexcept
    {
        std::vector<Scalar> res;
        for(size_t i = 0; i < this->N_m; i++){
            res.push_back(this->dr(i));
        }
        return res;
    }

    //! Get mesh dimensions
    /*!*************************************************************************
    * Returns the number of mesh points.
    ***************************************************************************/
    size_t dim() const noexcept { return N_m; }

    /*!*************************************************************************
    * Class for representing a generic point in a 1D mesh.
    ***************************************************************************/
    class mesh_point {
    private:
        const Mesh_base* m_m;
        Scalar i_m;
        mesh_point() = default;
    public:
        mesh_point(const Mesh_base* m, const Scalar i) : m_m(m), i_m(i) {}
        mesh_point(const mesh_point&) = default;
        mesh_point(mesh_point&&) = default;
        ~mesh_point() = default;

        mesh_point& operator=(const mesh_point&) = default;
        mesh_point& operator=(mesh_point&&) = default;

        Scalar r() const noexcept { return m_m->r(i_m);}
        Scalar r2() const noexcept { return m_m->r2(i_m);}
        Scalar dr() const noexcept { return m_m->dr(i_m);}
    };

    /*!*************************************************************************
    * Class for enabling iteration through a mesh, constant version
    ***************************************************************************/
    class const_iterator{
    private:
        const Mesh_base* m_m;
        Scalar i_m;
        const_iterator() = default;
    public:
        typedef const Scalar difference_type;
        typedef const mesh_point value_type;
        // typedef const mesh_point reference;
        // typedef const mesh_point pointer;
        // typedef std::random_access_iterator_tag iterator_category;

        const_iterator(const const_iterator&) = default;
        const_iterator(const_iterator&&) = default;
        const_iterator(const Mesh_base* m, const Scalar i) : m_m(m), i_m(i){}
        ~const_iterator() = default;
        const_iterator& operator=(const const_iterator&) = default;
        const_iterator& operator=(const_iterator&&) = default;

        value_type operator*(){return mesh_point(m_m, i_m);}

        bool operator==(const const_iterator& b)const noexcept{return i_m == b.i_m;}
        bool operator!=(const const_iterator& b)const noexcept{return !(*this == b);}
        bool operator>=(const const_iterator& b)const noexcept{return i_m >= b.i_m;}
        bool operator<=(const const_iterator& b)const noexcept{return i_m <= b.i_m;}
        bool operator>(const const_iterator& b)const noexcept{return !(*this <= b);}
        bool operator<(const const_iterator& b)const noexcept{return !(*this >= b);}

        const_iterator operator++(int) noexcept{auto res = *this; i_m++; return res;}
        const_iterator& operator++() noexcept{i_m++; return *this;}
        const_iterator operator--(int) noexcept{auto res = *this; i_m--; return res;}
        const_iterator& operator--() noexcept{i_m--; return *this;}

        const_iterator& operator+=(const difference_type n) noexcept{i_m += n; return *this;}
        const_iterator& operator-=(const difference_type n) noexcept{i_m -= n; return *this;}

        const_iterator operator+(const difference_type n)const noexcept{auto res = *this; return (res += n);}
        friend const_iterator operator+(const difference_type n, const_iterator it)
         noexcept { return it + n; }
        const_iterator operator-(const difference_type n)const noexcept{auto res = *this; return (res -= n);}
        difference_type operator-(const const_iterator& it)const noexcept{return this->i_m - it.i_m;}
    };

    /*!*************************************************************************
    * Class for enabling reverse iteration through a mesh, constant version
    ***************************************************************************/
    typedef std::reverse_iterator<const_iterator> const_reverse_iterator;

    const_iterator begin() const noexcept{return const_iterator(this, 0);}
    const_iterator cbegin() const noexcept{return const_iterator(this, 0);}

    const_iterator end() const noexcept{return const_iterator(this, static_cast<Scalar>(N_m));}
    const_iterator cend() const  noexcept{return const_iterator(this, static_cast<Scalar>(N_m));}

    const_reverse_iterator rbegin() const noexcept{return const_reverse_iterator(this->end());}
    const_reverse_iterator crbegin() const noexcept{return const_reverse_iterator(this->end());}

    const_reverse_iterator rend() const noexcept{return const_reverse_iterator(this->begin());}
    const_reverse_iterator crend() const noexcept{return const_reverse_iterator(this->cbegin());}
};

/*!*************************************************************************
* Class representing a D-dimensional linear mesh. Mesh points, \f$ \bar{r} \f$, are
* calculated as \f$ r_i = \alpha_i*x_i + \gamma_i, i = 0 \ldots D - 1 \f$
***************************************************************************/
template<size_t D, class Scalar = double>
class Linear_mesh : public Mesh_base<D, Scalar> {
private:
    Linear_mesh() = default;

public:
    using Arr = std::array<Scalar, D>;
    //! Constructor
    /*!*************************************************************************
    * Construct a linear mesh from two parameters, and the number of mesh points.
    * @param R_min Lower end points for the mesh.
    * @param R_max Upper end points for the mesh.
    * @param N Number of mesh points (including endpoints).
    ***************************************************************************/
    Linear_mesh(const Arr& R_min, const Arr& R_max, const std::array<size_t, D>& N)
     noexcept
    : Mesh_base<D, Scalar>(
        [this] (const Arr& r_min, const Arr& r_max,  const std::array<size_t, D>& n_steps)
        {
            Arr res;
            auto z3 = this->zip3(r_min,  r_max, n_steps);
            std::transform(std::begin(z3), std::end(z3), std::begin(res),
                [] (std::tuple<Scalar, Scalar, size_t> ti)
                {
#if __cplusplus >= 201703L
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

    using Mesh_base<D, Scalar>::r;
    using Mesh_base<D, Scalar>::r2;
    using Mesh_base<D, Scalar>::dr;



    //! Get position at coordinates
    /*!*************************************************************************
    * This function evaluates the mesh at the coordinates supplied, returning the
    * position, \f$ r_i(x_i) = \alpha_i*x_i + \gamma_i \f$.
    * @param x Array containing the coordinate to evaluate.
    * \return The position, r, of the point at index x.
    ***************************************************************************/
    Arr r(const Arr& x)
        const noexcept override
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
#if __cplusplus >= 201703L
            auto [xi, ai, bi, ci] = ti;
#else
            auto xi = std::get<0>(ti);
            auto ai = std::get<1>(ti);
            auto bi = std::get<2>(ti);
            auto ci = std::get<3>(ti);
#endif
            return ai*xi +0*bi + ci;
        });
#endif
        return res;
    }

    //! Get position squared at coordinates
    /*!*************************************************************************
    * This function evaluates the mesh at the coordinates supplied, returning an
    * Array containing the squared components of the position,
    * \f$ r2_i(x) = r_i(x)^2 \f$.
    * @param x Array containing the coordinate to evaluate.
    * \return The squared position, r2, of the point at index x.
    ***************************************************************************/
    Arr r2(const Arr& x)
        const noexcept override
    {
        Arr res, r = this->r(x);
        std::transform(std::begin(r), std::end(r), std::begin(res),
            [] (const Scalar s) {return s*s;});
        return res;
    }

    //! Get step size at coordinates
    /*!*************************************************************************
    * This function evaluates the mesh at the coordinates supplied, returning an
    * Array containing the step size at coordinate x,
    * \f$ dr(x)_i = \left(\left.\frac{dr}{dx}\right|_x\right)_i = \alpha_i \f$.
    * @param x Array containing the coordinate to evaluate.
    * \return The step size, \f$ \alpha_i \f$, of the point at index x.
    ***************************************************************************/
    Arr dr(const Arr& x)
        const noexcept override
    {
        return this->alpha_m;
    }
};

/*!*****************************************************************************
* Partial specialization Linear_mesh for one dimension,
* \f$ r = \alpha*x + \gamma \f$.
*******************************************************************************/
template<class Scalar>
class Linear_mesh<1, Scalar> : public Mesh_base<1, Scalar> {
private:
    Linear_mesh() = default;

public:
    //! Constructor
    /*!*************************************************************************
    * Construct a linear mesh from two parameters, and the number of mesh points.
    * @param R_min Lower end point for the mesh.
    * @param R_max Upper end point for the mesh.
    * @param N Number of mesh points (including endpoints).
    ***************************************************************************/
    Linear_mesh(const Scalar R_min, const Scalar R_max, const size_t N)
     noexcept
    : Mesh_base<1, Scalar>((R_max - R_min)/static_cast<double>(N - 1), 0, R_min, N)
    {}

    ~Linear_mesh() = default;

    Linear_mesh(const Linear_mesh&) = default;
    Linear_mesh(Linear_mesh&&) = default;

    Linear_mesh& operator=(const Linear_mesh&) = default;
    Linear_mesh& operator=(Linear_mesh&&) = default;

    using Mesh_base<1, Scalar>::r;
    using Mesh_base<1, Scalar>::r2;
    using Mesh_base<1, Scalar>::dr;

    //! Get position at coordinate
    /*!*************************************************************************
    * This function evaluates the mesh at the coordinate supplied, returning the
    * position, \f$ r(x) = \alpha*x + \gamma \f$.
    * @param x Coordinate to evaluate.
    * \return The position, r, of the point at index x.
    ***************************************************************************/
    Scalar r(const Scalar& x)
        const noexcept override
    {
        return x*this->alpha_m + this->gamma_m;
    }

    //! Get position squared at coordinate
    /*!*************************************************************************
    * This function evaluates the mesh at the coordinate squared supplied,
    * returning the squared position, \f$ r2(x) = (r(x))^2 \f$.
    * @param x Coordinate to evaluate.
    * \return The squared position, \f$ r2 \f$, of the point at index x.
    ***************************************************************************/
    Scalar r2(const Scalar& x)
        const noexcept override
    {
        return std::pow(this->r(x), 2);
    }

    //! Get step size at coordinate
    /*!*************************************************************************
    * This function evaluates the mesh at the coordinate squared supplied,
    * returning the step_size, \f$ dr(x) = \frac{dr}{dx} = \alpha \f$.
    * @param x Coordinate to evaluate.
    * \return The step size, \f$ \alpha \f$, at index x.
    ***************************************************************************/
    Scalar dr(const Scalar& x)
        const noexcept override
    {
        return this->alpha_m + 0*x;
    }
};


/*!*****************************************************************************
* Class representing a D-dimensional quadratic mesh. Mesh points,
* \f$ \bar{r} \f$, are calculated as \f$ r_i = \alpha_i*x_i^2 + \gamma_i \f$.
*******************************************************************************/
template<size_t D, class Scalar = double>
class Quadratic_mesh : public Mesh_base<D, Scalar> {
private:
    Quadratic_mesh() = default;
public:
    using Arr = std::array<Scalar, D>;
    //! Constructor
    /*!*************************************************************************
    * Construct a quadratic mesh from two parameters, and the number of mesh
    * points.
    * @param R_min Lower end points for the mesh.
    * @param R_max Upper end points for the mesh.
    * @param N Number of mesh points (including endpoints).
    ***************************************************************************/
    Quadratic_mesh(const Arr& R_min, const Arr& R_max, const std::array<size_t, D>& N)
     noexcept
    : Mesh_base<D, Scalar>(
        [this] (const Arr& r_min, const Arr& r_max, const std::array<size_t, D>& n_steps)
        {
            Arr res;
            auto z3 = this->zip3(r_min, r_max, n_steps);
            std::transform(std::begin(z3), std::end(z3), std::begin(res),
                [] (const std::tuple<Scalar, Scalar, size_t>& ti)
                {
#if __cplusplus >= 201703L
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

    //! Get position at coordinates
    /*!*************************************************************************
    * This function evaluates the mesh at the coordinates supplied, returning the
    * position, \f$ r_i(x) = \alpha_i*x_i*x_i + \gamma_i \f$.
    * @param x Array containing the coordinate to evaluate.
    * \return The position, r, of the point at index x.
    ***************************************************************************/
    Arr r(const Arr& x)
        const noexcept override
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
#if __cplusplus >= 201703L
                auto [xi, a, b, c] = ti;
#else
                auto xi = std::get<0>(ti);
                auto a = std::get<1>(ti);
                auto b = std::get<2>(ti);
                auto c = std::get<3>(ti);
#endif
                return sgn(xi)*a*std::pow(xi, 2) + c;
            });
#endif
        return res;
    }

    //! Get position squared at coordinate
    /*!*************************************************************************
    * This function evaluates the mesh at the coordinate squared supplied,
    * returning the squared position, \f$ r2(x) = (r(x))^2 \f$.
    * @param x Coordinate to evaluate.
    * \return The squared position, \f$ r2 \f$, of the point at index x.
    ***************************************************************************/
    Arr r2(const Arr& x)
        const noexcept override
    {
        Arr res, r = this->r(x);
        std::transform(std::begin(r), std::end(r), std::begin(res), [] (const Scalar& s) {return s*s;});
        return res;
    }

    //! Get step size at coordinates
    /*!*************************************************************************
    * This function evaluates the mesh at the coordinates supplied, returning an
    * Array containing the step size at coordinate x,
    * \f$ dr(x_i)_i = \left(\left.\frac{dr}{dx}\right|_x\right)_i = 2x_i\alpha_i \f$.
    * @param x Array containing the coordinate to evaluate.
    * \return The step size, \f$ 2x_i\alpha_i \f$, of the point at index x.
    ***************************************************************************/
    Arr dr(const Arr& x)
        const noexcept override
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
#if __cplusplus >= 201703L
                auto [xi, a, b, c] = ti;
#else
                auto xi = std::get<0>(ti);
                auto a = std::get<1>(ti);
                auto b = std::get<2>(ti);
                auto c = std::get<3>(ti);
#endif
                return sgn(xi)*2*a*xi + b*0 + c*0;
            });
#endif
        return res;
    }
};

/*!*****************************************************************************
* Quadratic mesh in one dimension, \f$ r = \alpha*x^2 + \gamma \f$.
*******************************************************************************/
template<class Scalar>
class Quadratic_mesh<1, Scalar> : public Mesh_base<1, Scalar> {
private:
    Quadratic_mesh() = default;
public:
    //! Constructor
    /*!*************************************************************************
    * Construct a 1D-quadratic mesh from two parameters, and the number of mesh
    * points.
    * @param R_min Lower end point for the mesh.
    * @param R_max Upper end point for the mesh.
    * @param N Number of mesh points (including endpoints).
    ***************************************************************************/
    Quadratic_mesh(const Scalar R_min, const Scalar R_max, const size_t N)
     noexcept
    : Mesh_base<1, Scalar>((R_max - R_min)/std::pow(N - 1, 2), 0, R_min, N )
    {}

    ~Quadratic_mesh() = default;

    Quadratic_mesh(const Quadratic_mesh&) = default;
    Quadratic_mesh(Quadratic_mesh&&) = default;

    Quadratic_mesh& operator=(const Quadratic_mesh&) = default;
    Quadratic_mesh& operator=(Quadratic_mesh&&) = default;

    using Mesh_base<1, Scalar>::r;
    using Mesh_base<1, Scalar>::r2;
    using Mesh_base<1, Scalar>::dr;

    //! Get position at coordinates
    /*!*************************************************************************
    * This function evaluates the mesh at the coordinates supplied, returning the
    * position, \f$ r(x) = \alpha*x^2 + \gamma \f$.
    * @param x Coordinate to evaluate.
    * \return The position, r, of the point at index x.
    ***************************************************************************/
    Scalar r(const Scalar& x)
        const noexcept override
    {
        return sgn(x)*std::pow(x, 2)*this->alpha_m + this->gamma_m;
    }

    //! Get position squared at coordinate
    /*!*************************************************************************
    * This function evaluates the mesh at the coordinate squared supplied,
    * returning the squared position, \f$ r2(x) = (r(x))^2 \f$.
    * @param x Coordinate to evaluate.
    * \return The squared position, \f$ r2 \f$, of the point at index x.
    ***************************************************************************/
    Scalar r2(const Scalar& x)
        const noexcept override
    {
        return std::pow(this->r(x), 2);
    }

    //! Get step size at coordinates
    /*!*************************************************************************
    * This function evaluates the mesh at the coordinates supplied, returning
    * the step size at coordinate x,
    * \f$ dr(x) = \left.\frac{dr}{dx}\right|_x\ = 2x\alpha \f$.
    * @param x Array containing the coordinate to evaluate.
    * \return The step size, \f$ 2x\alpha \f$, of the point at index x.
    ***************************************************************************/
    Scalar dr(const Scalar& x)
        const noexcept override
    {
        return 2*this->alpha_m*sgn(x)*x;
    }
};

/*!*****************************************************************************
* Exponential mesh in arbitrary dimensions,
* \f$ r_i = \alpha_i(e^{\beta_i x} - 1) + \gamma_i \f$.
*******************************************************************************/
template<size_t D, class Scalar = double>
class Exponential_mesh : public Mesh_base<D, Scalar> {
private:
    Exponential_mesh() = default;
public:
    using Arr = std::array<Scalar, D>;
    //! Constructor
    /*!*************************************************************************
    * Construct a exponential mesh from three parameters, and the number of mesh
    * points.
    * @param R_min Lower end points for the mesh.
    * @param R_max Upper end points for the mesh.
    * @param beta Array containing the \f$ \beta_i \f$ for the mesh.
    * @param N Number of mesh points (including endpoints).
    ***************************************************************************/
    Exponential_mesh(const Arr& R_min, const Arr& R_max, const Arr& beta, const std::array<size_t, D>& N)
     noexcept
    : Mesh_base<D, Scalar>(
        [this] (const Arr& r_min, const Arr& r_max, const Arr& beta_val, const std::array<size_t, D>& n_steps)
        {
            Arr res;
            auto z4 = this->zip4(r_min, r_max, beta_val, n_steps);
            std::transform(std::begin(z4), std::end(z4), std::begin(res),
                [] (const std::tuple<Scalar, Scalar, Scalar, size_t> ti)
                {
#if __cplusplus >= 201703L
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

    //! Get position at coordinates
    /*!*************************************************************************
    * This function evaluates the mesh at the coordinates supplied, returning the
    * position, \f$ r_i(x_i) = \alpha_i(e^{\beta_i x_i} - 1) + \gamma_i \f$.
    * @param x Array containing the coordinate to evaluate.
    * \return The position, r, of the point at index x.
    ***************************************************************************/
    Arr r(const Arr& x)
        const noexcept override
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
#if __cplusplus >= 201703L
                auto [xi, a, b, c] = ti;
#else
                auto xi = std::get<0>(ti);
                auto a = std::get<1>(ti);
                auto b = std::get<2>(ti);
                auto c = std::get<3>(ti);
#endif
                return sgn(xi)*a*(std::exp(sgn(xi)*xi*b) - 1) + c;
            });
#endif
        return res;
    }

    //! Get position squared at coordinate
    /*!*************************************************************************
    * This function evaluates the mesh at the coordinate squared supplied,
    * returning the squared position, \f$ r2_i(x_i) = (r(x_i))^2 \f$.
    * @param x Coordinate to evaluate.
    * \return The squared position, \f$ r2 \f$, of the point at index x.
    ***************************************************************************/
    Arr r2(const Arr& x)
        const noexcept override
    {
        Arr res, r = this->r(x);
        std::transform(std::begin(r), std::end(r), std::begin(res), [] (const Scalar& s) {return s*s;});
        return res;
    }

    //! Get step size at coordinates
    /*!*************************************************************************
    * This function evaluates the mesh at the coordinates supplied, returning
    * the step size at coordinate x,
    * \f$ dr_i(x_i) = \left(\left.\frac{dr}{dx}\right|_x\right)_i = \alpha_i\beta_i e^{\beta_i x_i} \f$.
    * @param x Array containing the coordinate to evaluate.
    * \return The step size, \f$ \alpha_i\beta_i e^{\beta_i x_i} \f$, of the
    * point at index x.
    ***************************************************************************/
    Arr dr(const Arr& x)
        const noexcept override
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
#if __cplusplus >= 201703L
                auto [ri, a, b, c] = ti;
#else
                auto ri = std::get<0>(ti);
                auto a = std::get<1>(ti);
                auto b = std::get<2>(ti);
                auto c = std::get<3>(ti);
#endif
                return b*(ri + a - c);
            });
#endif
        return res;
    }
};

/*!*****************************************************************************
* Exponential mesh in one dimension, \f$ r = \alpha(e^{\beta x} - 1) + \gamma \f$.
*******************************************************************************/
template<class Scalar>
class Exponential_mesh<1, Scalar> : public Mesh_base<1, Scalar> {
private:
    Exponential_mesh() = default;
public:
    //! Constructor
    /*!*************************************************************************
    * Construct a exponential mesh from three parameters, and the number of mesh
    * points.
    * @param R_min Lower end point for the mesh.
    * @param R_max Upper end point for the mesh.
    * @param beta The parameter \f$ \beta \f$ for the mesh.
    * @param N Number of mesh points (including endpoints).
    ***************************************************************************/
    Exponential_mesh(const Scalar R_min, const Scalar R_max, const Scalar beta, const size_t N)
     noexcept
    : Mesh_base<1, Scalar>((R_max - R_min)/std::expm1(beta*static_cast<double>(N-1)), beta, R_min, N)
    {}

    ~Exponential_mesh() = default;

    Exponential_mesh(const Exponential_mesh&) = default;
    Exponential_mesh(Exponential_mesh&&) = default;

    Exponential_mesh& operator=(const Exponential_mesh&) = default;
    Exponential_mesh& operator=(Exponential_mesh&&) = default;

    using Mesh_base<1, Scalar>::r;
    using Mesh_base<1, Scalar>::r2;
    using Mesh_base<1, Scalar>::dr;

    //! Get position at coordinates
    /*!*************************************************************************
    * This function evaluates the mesh at the coordinates supplied, returning the
    * position, \f$ r(x) = \alpha(e^{\beta x} - 1) + \gamma \f$.
    * @param x Array containing the coordinate to evaluate.
    * \return The position, r, of the point at index x.
    ***************************************************************************/
    Scalar r(const Scalar& x)
        const noexcept override
    {
        Scalar a = this->alpha_m, b = this->beta_m, c = this->gamma_m;
        return sgn(x)*a*std::expm1(sgn(x)*b*x) + c;
    }

    //! Get position squared at coordinate
    /*!*************************************************************************
    * This function evaluates the mesh at the coordinate squared supplied,
    * returning the squared position, \f$ r2(x) = (r(x))^2 \f$.
    * @param x Coordinate to evaluate.
    * \return The squared position, \f$ r2 \f$, of the point at index x.
    ***************************************************************************/
    Scalar r2(const Scalar& x)
        const noexcept override
    {
        return std::pow(this->r(sgn(x)*x), 2);
    }

    //! Get step size at coordinates
    /*!*************************************************************************
    * This function evaluates the mesh at the coordinates supplied, returning
    * the step size at coordinate x,
    * \f$ dr(x) = \left.\frac{dr}{dx}\right|_x = \alpha\beta e^{\beta x} \f$.
    * @param x Array containing the coordinate to evaluate.
    * \return The step size, \f$ \alpha\beta e^{\beta x} \f$, of the point at
    * index x.
    ***************************************************************************/
    Scalar dr(const Scalar& x)
        const noexcept override
    {
        Scalar a = this->alpha_m, b = this->beta_m, c = this->gamma_m;
        return b*(this->r(sgn(x)*x) + a - c);
    }

    Scalar A()
        const noexcept
    {
        return this->alpha_m;
    }

    Scalar B()
        const noexcept
    {
        return this->beta_m;
    }

};

/*!*****************************************************************************
* @}
*******************************************************************************/

#endif //NUMERICAL_MESH_LIB_H
