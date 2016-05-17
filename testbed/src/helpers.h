// Helper functions for BD/MD

#ifndef BUFFMD_HELPERS_H_
#define BUFFMD_HELPERS_H_

#include <tuple>

// Template functions for everybody to use for tuple stuff
// Taken from Stack Overflow entirely
// XXX: Not C++14 compatible!!!!!!
namespace detail {
    template<bool b, typename T=void>
    using EnableIf = typename std::enable_if<b,T>::type;
    
    template<int Index, template<typename T>class Search, int Which, typename, class First, class... Types>
    struct get_internal:
        get_internal<Index + 1, Search, Which, void, Types...>
    {};

    template<int Index, template<typename T>class Search, int Which, class First, class... Types>
    struct get_internal<Index, Search, Which, EnableIf<!Search<First>::value>, First, Types...>:
        get_internal<Index + 1, Search, Which, void, Types...>
    {};
    template<int Index, template<typename T>class Search, int Which, class First, class... Types>
    struct get_internal<Index, Search, Which, EnableIf<Search<First>::value>, First, Types...>:
        get_internal<Index + 1, Search, Which-1, void, Types...>
    {};
    template<int Index, template<typename T>class Search, class First, class... Types>
    struct get_internal<Index, Search, 0, EnableIf<Search<First>::value>, First, Types...>:
        std::integral_constant<int, Index>
    {};

    template<template<typename>class Test, int Which=0, class... Types>
    auto get(std::tuple<Types...>& tuple)->
    decltype(std::get<get_internal<0,Test,Which,void,Types...>::value>(tuple))
    {
        return std::get<get_internal<0,Test,Which,void,Types...>::value>(tuple);
    }
    template<template<typename>class Test, int Which=0, class... Types>
    auto get(std::tuple<Types...> const& tuple)->
    decltype(std::get<get_internal<0,Test,Which,void,Types...>::value>(tuple))
    {
        return std::get<get_internal<0,Test,Which,void,Types...>::value>(tuple);
    }
    template<template<typename>class Test, int Which=0, class... Types>
    auto get(std::tuple<Types...>&& tuple)->
    decltype(std::move(std::get<get_internal<0,Test,Which,void,Types...>::value>(tuple)))
    {
        return std::move(std::get<get_internal<0,Test,Which,void,Types...>::value>(tuple));
    }

    template<typename T>
    struct is_type {
    template<typename U>
    using test = std::is_same<T,U>;
    };

    template<class T, int Which=0, class... Types>
    T& get(std::tuple<Types...>& tuple)
    {
        return get<is_type<T>::template test,Which>(tuple);
    }
    template<class T, int Which=0, class... Types>
    T const& get(std::tuple<Types...> const& tuple)
    {
        return get<is_type<T>::template test,Which>(tuple);
    }
    template<class T, int Which=0, class... Types>
    T&& get(std::tuple<Types...>&& tuple)
    {
        return std::move(get<is_type<T>::template test,Which>(tuple));
    }
}

namespace buffmd {
    __attribute__((always_inline,pure))
    static inline float pbc_float(float x, const float boxby2, const float box) {
        while (x >  boxby2) x -= box;
        while (x < -boxby2) x += box;
        return x;
    }

    __attribute__((always_inline,pure))
    static inline double pbc(double x, const double boxby2, const double box) {
        while (x >  boxby2) x -= box;
        while (x < -boxby2) x += box;
        return x;
    }
    
    __attribute__((always_inline,pure))
    static inline double pbc2(double x, const double box) {
        while (x > box) x -= box;
        while (x < 0  ) x += box;
        return x;
    }

    __attribute((always_inline))
    static inline void azzero(double *d, const int n) {
        for (int i = 0; i < n; ++i) {
            d[i] = 0.0;
        }
    }

    template<class... T> void unused(T&&...) {}
    
    // From bithacks online
    __attribute__((always_inline))
    static inline unsigned int nextpow2(unsigned int v) {
        v--;
        
        v |= v >> 1;
        v |= v >> 2;
        v |= v >> 4;
        v |= v >> 8;
        v |= v >> 16;
        v++;
        
        return v;
    }
    
    __attribute__((always_inline))
    static inline int cell_vec_to_linear(int cx, int cy, int cz, int nc[3]) {
        return cx + cy*nc[0] + cz*nc[0]*nc[1];
    }
    // cidx = buffmd::cell_vec_to_linear((int[]){cx,cy,cz}, T_);
    
    __attribute__((always_inline))
    static inline void cell_linear_to_vec(int cidx, int nc[3], int* cx) {
        cx[0] = cidx % nc[0];
        cx[1] = (cidx / nc[0]) % nc[1];
        cx[2] = cidx / (nc[0] * nc[1]);
    }
    // buffmd::cell_linear_to_vec(cidx, T_, c_test);
}

#endif
