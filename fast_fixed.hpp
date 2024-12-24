#pragma once

#include <iostream>
#include <type_traits>
#include <cstdint>
#include <cmath>

template<size_t N, size_t K>
struct Fixed;

template<size_t N, size_t K>
struct FastFixed {
    static_assert(N > K, "N must be greater than K");

    static constexpr size_t bits = N;
    static constexpr size_t frac = K;

    using StorageType = typename std::conditional<(N <= 8), int_fast8_t,
        typename std::conditional<(N <= 16), int_fast16_t,
            typename std::conditional<(N <= 32), int_fast32_t, int_fast64_t>::type
        >::type
    >::type;

    StorageType v;

    constexpr FastFixed(int v = 0): v(static_cast<StorageType>(v) << K) {}
    constexpr FastFixed(float f): v(f * (StorageType(1) << K)) {}
    constexpr FastFixed(double f): v(f * (StorageType(1) << K)) {}

    static constexpr FastFixed from_raw(StorageType x) {
        FastFixed ret;
        ret.v = x;
        return ret;
    }

    auto operator<=>(const FastFixed&) const = default;
    bool operator==(const FastFixed&) const = default;

    explicit operator float() const { return v / float(StorageType(1) << K); }
    explicit operator double() const { return v / double(StorageType(1) << K); }

    friend FastFixed operator/(FastFixed a, int b) {
        return FastFixed::from_raw(a.v / b);
    }

    friend FastFixed operator*(FastFixed a, int b) {
        return FastFixed::from_raw(a.v * b);
    }

    friend FastFixed operator*(int a, FastFixed b) {
        return b * a;
    }

    template<size_t N2, size_t K2>
    explicit operator FastFixed<N2,K2>() const {
        if constexpr (K2 >= K) {
            return FastFixed<N2,K2>::from_raw(static_cast<typename FastFixed<N2,K2>::StorageType>(v) << (K2 - K));
        } else {
            constexpr size_t shift = K - K2;
            if constexpr (shift >= N2) {
                auto temp = v >> (shift - N2 + 1);
                return FastFixed<N2,K2>::from_raw(static_cast<typename FastFixed<N2,K2>::StorageType>(temp) >> 1);
            } else {
                return FastFixed<N2,K2>::from_raw(static_cast<typename FastFixed<N2,K2>::StorageType>(v) >> shift);
            }
        }
    }

    template<size_t N2, size_t K2>
    explicit operator Fixed<N2,K2>() const {
        if constexpr (K2 >= K) {
            return Fixed<N2,K2>::from_raw(static_cast<typename Fixed<N2,K2>::StorageType>(v) << (K2 - K));
        } else {
            constexpr size_t shift = K - K2;
            if constexpr (shift >= N2) {
                auto temp = v >> (shift - N2 + 1);
                return Fixed<N2,K2>::from_raw(static_cast<typename Fixed<N2,K2>::StorageType>(temp) >> 1);
            } else {
                return Fixed<N2,K2>::from_raw(static_cast<typename Fixed<N2,K2>::StorageType>(v) >> shift);
            }
        }
    }
    public:
        static constexpr std::size_t kNValue = N;
        static constexpr std::size_t kKValue = K;
        static constexpr bool kFast          = true;

};

template<size_t N, size_t K>
FastFixed<N,K> operator+(FastFixed<N,K> a, FastFixed<N,K> b) {
    return FastFixed<N,K>::from_raw(a.v + b.v);
}

template<size_t N, size_t K>
FastFixed<N,K> operator-(FastFixed<N,K> a, FastFixed<N,K> b) {
    return FastFixed<N,K>::from_raw(a.v - b.v);
}

template<size_t N, size_t K>
FastFixed<N,K> operator*(FastFixed<N,K> a, FastFixed<N,K> b) {
    using ST = typename FastFixed<N,K>::StorageType;
    if constexpr (N <= 32) {
        return FastFixed<N,K>::from_raw((static_cast<int_fast64_t>(a.v) * b.v) >> K);
    } else {
        ST high = (a.v >> K) * b.v;
        ST low = (a.v & ((ST(1) << K) - 1)) * b.v >> K;
        return FastFixed<N,K>::from_raw(high + low);
    }
}

template<size_t N, size_t K>
FastFixed<N,K> operator/(FastFixed<N,K> a, FastFixed<N,K> b) {
    using ST = typename FastFixed<N,K>::StorageType;
    if constexpr (N <= 32) {
        return FastFixed<N,K>::from_raw((static_cast<int_fast64_t>(a.v) << K) / b.v);
    } else {
        return FastFixed<N,K>::from_raw((a.v << K) / b.v);
    }
}

template<size_t N, size_t K>
FastFixed<N,K> &operator+=(FastFixed<N,K> &a, FastFixed<N,K> b) { return a = a + b; }

template<size_t N, size_t K>
FastFixed<N,K> &operator-=(FastFixed<N,K> &a, FastFixed<N,K> b) { return a = a - b; }

template<size_t N, size_t K>
FastFixed<N,K> &operator*=(FastFixed<N,K> &a, FastFixed<N,K> b) { return a = a * b; }

template<size_t N, size_t K>
FastFixed<N,K> &operator/=(FastFixed<N,K> &a, FastFixed<N,K> b) { return a = a / b; }

template<size_t N, size_t K>
FastFixed<N,K> operator-(FastFixed<N,K> x) { return FastFixed<N,K>::from_raw(-x.v); }

template<size_t N, size_t K>
FastFixed<N,K> abs(FastFixed<N,K> x) {
    FastFixed<N,K> ret = x;
    if (ret.v < 0) ret.v = -ret.v;
    return ret;
}

template<size_t N, size_t K>
std::ostream &operator<<(std::ostream &out, FastFixed<N,K> x) {
    return out << static_cast<double>(x);
}
