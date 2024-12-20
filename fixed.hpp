// fixed.hpp

#ifndef FIXED_HPP
#define FIXED_HPP

#include <iostream>
#include <type_traits>
#include <cstdint>
#include <cassert>

// Helper to select the appropriate integer type based on bit size
template<size_t N>
struct IntTypeSelector;

template<>
struct IntTypeSelector<8> {
    using type = int8_t;
};

template<>
struct IntTypeSelector<16> {
    using type = int16_t;
};

template<>
struct IntTypeSelector<32> {
    using type = int32_t;
};

template<>
struct IntTypeSelector<64> {
    using type = int64_t;
};

template<size_t N, size_t K>
class Fixed {
    static_assert(K < N, "Fractional bits must be less than total bits");
    static_assert(N == 8 || N == 16 || N == 32 || N == 64, "N must be 8, 16, 32, or 64");

    using StorageType = typename IntTypeSelector<N>::type;
    StorageType value;

public:
    constexpr Fixed(float f) : value(static_cast<StorageType>(f * (1 << K))) {}
    constexpr Fixed(double d) : value(static_cast<StorageType>(d * (1 << K))) {}

    static constexpr Fixed from_raw(StorageType raw) {
        Fixed fp;
        fp.value = raw;
        return fp;
    }

    friend std::ostream& operator<<(std::ostream& os, const Fixed& fp) {
        return os << static_cast<double>(fp.value) / (1 << K);
    }

    // Arithmetic operations
    friend Fixed operator+(Fixed a, Fixed b) {
        return Fixed::from_raw(a.value + b.value);
    }

    friend Fixed operator-(Fixed a, Fixed b) {
        return Fixed::from_raw(a.value - b.value);
    }

    friend Fixed operator*(Fixed a, Fixed b) {
        return Fixed::from_raw((static_cast<int64_t>(a.value) * b.value) >> K);
    }

    friend Fixed operator/(Fixed a, Fixed b) {
        return Fixed::from_raw((static_cast<int64_t>(a.value) << K) / b.value);
    }

    Fixed& operator+=(Fixed b) {
        value += b.value;
        return *this;
    }

    Fixed& operator-=(Fixed b) {
        value -= b.value;
        return *this;
    }

    Fixed& operator*=(Fixed b) {
        value = (static_cast<int64_t>(value) * b.value) >> K;
        return *this;
    }

    Fixed& operator/=(Fixed b) {
        value = (static_cast<int64_t>(value) << K) / b.value;
        return *this;
    }
};

#endif // FIXED_HPP