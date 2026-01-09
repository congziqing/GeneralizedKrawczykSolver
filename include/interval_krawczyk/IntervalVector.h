#pragma once

#include <array>
#include <cmath>
#include <iostream>
#include <algorithm>
#include <stdexcept>
#include <Eigen/Dense>
#include "KaucherInterval.h"

namespace ik {

// Fixed-size interval vector class
template<size_t N>
class IntervalVector
{
public:
    std::array<KaucherInterval, N> data;

    static constexpr size_t Size = N;

    IntervalVector()
    {
        data.fill(KaucherInterval());
    }

    IntervalVector(std::initializer_list<KaucherInterval> init)
    {
        size_t i = 0;
        for (auto& val : init)
        {
            if (i < N) data[i++] = val;
        }
    }

    static IntervalVector<N> zeros()
    {
        IntervalVector<N> result;
        result.data.fill(KaucherInterval(0.0, 0.0));
        return result;
    }

    static IntervalVector<N> ones()
    {
        IntervalVector<N> result;
        result.data.fill(KaucherInterval(1.0, 1.0));
        return result;
    }

    static IntervalVector<N> constants(double c)
    {
        IntervalVector<N> result;
        result.data.fill(KaucherInterval(c, c));
        return result;
    }

    static IntervalVector<N> constants(double lower, double upper)
    {
        IntervalVector<N> result;
        result.data.fill(KaucherInterval(lower, upper));
        return result;
    }

    static IntervalVector<N> fromPoint(const Eigen::VectorXd& v)
    {
        IntervalVector<N> result;
        for (size_t i = 0; i < N; ++i)
        {
            result[i] = KaucherInterval(v(static_cast<int>(i)), v(static_cast<int>(i)));
        }
        return result;
    }

    static IntervalVector<N> fromBounds(const Eigen::VectorXd& lower, const Eigen::VectorXd& upper)
    {
        IntervalVector<N> result;
        for (size_t i = 0; i < N; ++i)
        {
            result[i] = KaucherInterval(lower(static_cast<int>(i)), upper(static_cast<int>(i)));
        }
        return result;
    }

    KaucherInterval& operator[](size_t i)
    {
        return data[i];
    }

    const KaucherInterval& operator[](size_t i) const
    {
        return data[i];
    }

    IntervalVector<N> operator+(const IntervalVector<N>& other) const
    {
        IntervalVector<N> result;
        for (size_t i = 0; i < N; ++i)
            result[i] = data[i] + other[i];
        return result;
    }

    IntervalVector<N> operator-(const IntervalVector<N>& other) const
    {
        IntervalVector<N> result;
        for (size_t i = 0; i < N; ++i)
            result[i] = data[i] - other[i];
        return result;
    }

    IntervalVector<N> operator-() const
    {
        IntervalVector<N> result;
        for (size_t i = 0; i < N; ++i)
            result[i] = -data[i];
        return result;
    }

    IntervalVector<N> operator*(double s) const
    {
        IntervalVector<N> result;
        for (size_t i = 0; i < N; ++i)
            result[i] = data[i] * s;
        return result;
    }

    friend IntervalVector<N> operator*(double s, const IntervalVector<N>& v)
    {
        return v * s;
    }

    IntervalVector<N> operator/(double s) const
    {
        if (s == 0.0)
            throw std::runtime_error("Division by zero");
        return (*this) * (1.0 / s);
    }

    IntervalVector<N>& operator+=(const IntervalVector<N>& other)
    {
        for (size_t i = 0; i < N; ++i)
            data[i] += other[i];
        return *this;
    }

    IntervalVector<N>& operator-=(const IntervalVector<N>& other)
    {
        for (size_t i = 0; i < N; ++i)
            data[i] -= other[i];
        return *this;
    }

    IntervalVector<N>& operator*=(double s)
    {
        for (size_t i = 0; i < N; ++i)
            data[i] *= s;
        return *this;
    }

    IntervalVector<N> elementwiseProduct(const IntervalVector<N>& other) const
    {
        IntervalVector<N> result;
        for (size_t i = 0; i < N; ++i)
            result[i] = data[i] * other[i];
        return result;
    }

    IntervalVector<N> elementwiseDivision(const IntervalVector<N>& other) const
    {
        IntervalVector<N> result;
        for (size_t i = 0; i < N; ++i)
            result[i] = data[i] / other[i];
        return result;
    }

    IntervalVector<N> intersection(const IntervalVector<N>& other) const
    {
        IntervalVector<N> result;
        for (size_t i = 0; i < N; ++i)
            result[i] = data[i].meet(other[i]);
        return result;
    }

    static IntervalVector<N> hull(const IntervalVector<N>& a, const IntervalVector<N>& b)
    {
        IntervalVector<N> result;
        for (size_t i = 0; i < N; ++i)
            result[i] = a[i].join(b[i]);
        return result;
    }

    bool contains(const IntervalVector<N>& other) const
    {
        for (size_t i = 0; i < N; ++i)
            if (!data[i].contains(other[i]))
                return false;
        return true;
    }

    bool contains(const Eigen::VectorXd& point) const
    {
        for (size_t i = 0; i < N; ++i)
        {
            if (!data[i].contains(point(static_cast<int>(i))))
                return false;
        }
        return true;
    }

    bool empty() const
    {
        return data[0].isEmpty();
    }

    double width() const
    {
        double max_width = 0.0;
        for (size_t i = 0; i < N; ++i)
        {
            double w = data[i].magnitude();
            if (w > max_width)
                max_width = w;
        }
        return max_width;
    }

    double maxWidth() const
    {
        return width();
    }

    IntervalVector<N> middle() const
    {
        IntervalVector<N> result;
        for (size_t i = 0; i < N; ++i)
            result[i] = KaucherInterval(data[i].middle());
        return result;
    }

    Eigen::VectorXd midpoint() const
    {
        Eigen::VectorXd result(N);
        for (size_t i = 0; i < N; ++i)
            result(static_cast<int>(i)) = data[i].middle();
        return result;
    }

    Eigen::VectorXd lower() const
    {
        Eigen::VectorXd result(N);
        for (size_t i = 0; i < N; ++i)
            result(static_cast<int>(i)) = data[i].lower();
        return result;
    }

    Eigen::VectorXd upper() const
    {
        Eigen::VectorXd result(N);
        for (size_t i = 0; i < N; ++i)
            result(static_cast<int>(i)) = data[i].upper();
        return result;
    }

    void expandToInclude(const IntervalVector<N>& other)
    {
        for (size_t i = 0; i < N; ++i)
            data[i].expandToInclude(other[i]);
    }

    IntervalVector<N> copy() const
    {
        IntervalVector<N> result;
        for (size_t i = 0; i < N; ++i)
            result[i] = data[i];
        return result;
    }

    void print(std::ostream& os = std::cout, const std::string& delimiter = ", ") const
    {
        os << "[";
        for (size_t i = 0; i < N; ++i)
        {
            os << data[i];
            if (i < N - 1)
                os << delimiter;
        }
        os << "]";
    }

    friend std::ostream& operator<<(std::ostream& os, const IntervalVector<N>& v)
    {
        v.print(os);
        return os;
    }

    Eigen::VectorXd toEigen() const
    {
        return midpoint();
    }
};

namespace interval
{
    template<size_t N>
    inline IntervalVector<N> abs(const IntervalVector<N>& v)
    {
        IntervalVector<N> result;
        for (size_t i = 0; i < N; ++i)
            result[i] = KaucherInterval(v[i].absoluteValue());
        return result;
    }

    template<size_t N>
    inline IntervalVector<N> sqrt(const IntervalVector<N>& v)
    {
        IntervalVector<N> result;
        for (size_t i = 0; i < N; ++i)
            result[i] = interval::sqrt(v[i]);
        return result;
    }

    template<size_t N>
    inline IntervalVector<N> exp(const IntervalVector<N>& v)
    {
        IntervalVector<N> result;
        for (size_t i = 0; i < N; ++i)
            result[i] = interval::exp(v[i]);
        return result;
    }

    template<size_t N>
    inline IntervalVector<N> log(const IntervalVector<N>& v)
    {
        IntervalVector<N> result;
        for (size_t i = 0; i < N; ++i)
            result[i] = interval::log(v[i]);
        return result;
    }

    template<size_t N>
    inline IntervalVector<N> sin(const IntervalVector<N>& v)
    {
        IntervalVector<N> result;
        for (size_t i = 0; i < N; ++i)
            result[i] = interval::sin(v[i]);
        return result;
    }

    template<size_t N>
    inline IntervalVector<N> cos(const IntervalVector<N>& v)
    {
        IntervalVector<N> result;
        for (size_t i = 0; i < N; ++i)
            result[i] = interval::cos(v[i]);
        return result;
    }

    template<size_t N>
    inline IntervalVector<N> pow(const IntervalVector<N>& v, double p)
    {
        IntervalVector<N> result;
        for (size_t i = 0; i < N; ++i)
            result[i] = interval::pow(v[i], p);
        return result;
    }
}
} // namespace ik