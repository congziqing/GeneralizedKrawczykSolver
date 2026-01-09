#pragma once

#define _USE_MATH_DEFINES

#include <algorithm>
#include <cmath>
#include <iostream>
#include <limits>
#include <vector>
#include <stdexcept>

namespace ik {

class KaucherInterval
{
public:
    double start;
    double end;

    KaucherInterval()
        : start(std::numeric_limits<double>::infinity()), end(-std::numeric_limits<double>::infinity())
    {
    }

    KaucherInterval(double s, double e) : start(s), end(e)
    {
    }

    explicit KaucherInterval(double x) : start(x), end(x)
    {
    }

    bool isProper() const
    {
        return start <= end;
    }

    bool isProper(const double tol) const
    {
        return start <= (end + tol);
    }

    bool isImproper() const
    {
        return start > end;
    }

    bool isImproper(const double tol) const
    {
        return start > (end + tol);
    }

    bool isPoint() const
    {
        return start == end;
    }

    bool isInfinite() const
    {
        return std::isinf(start) || std::isinf(end);
    }

    bool isEmpty() const
    {
        return start > end && std::isinf(start);
    }

    double width() const
    {
        return end - start;
    }

    double magnitude() const
    {
        return std::abs(end - start);
    }

    double absoluteValue() const
    {
        if (isEmpty()) return 0.0;
        if (isImproper()) return std::max(std::abs(start), std::abs(end));
        if (start >= 0) return end;
        if (end <= 0) return -start;
        return std::max(std::abs(start), std::abs(end));
    }

    double middle() const
    {
        return 0.5 * (start + end);
    }

    double lower() const
    {
        return std::min(start, end);
    }

    double upper() const
    {
        return std::max(start, end);
    }

    KaucherInterval dual() const
    {
        return KaucherInterval(end, start);
    }

    KaucherInterval proper() const
    {
        if (isProper())
            return *this;
        return dual();
    }

    void makeProper()
    {
        if (start > end)
        {
            std::swap(start, end);
        }
    }

    KaucherInterval meet(const KaucherInterval& other) const
    {
        if (isInfinite())
            return other;
        if (other.isInfinite())
            return *this;

        return KaucherInterval(std::max(start, other.start), std::min(end, other.end));
    }

    KaucherInterval join(const KaucherInterval& other) const
    {
        if (start > end && std::isinf(start))
            return other;
        if (other.start > other.end && std::isinf(other.start))
            return *this;

        return KaucherInterval(std::min(start, other.start), std::max(end, other.end));
    }

    void expandToInclude(double x)
    {
        if (start > end && std::isinf(start))
        {
            start = end = x;
        }
        else
        {
            start = std::min(start, x);
            end = std::max(end, x);
        }
    }

    void expandToInclude(const KaucherInterval& other)
    {
        if (isEmpty())
        {
            start = other.start;
            end = other.end;
        }
        else if (!other.isEmpty())
        {
            start = std::min(start, other.start);
            end = std::max(end, other.end);
        }
    }

    bool contains(double x) const
    {
        if (isProper())
            return x >= start && x <= end;
        return x >= end && x <= start;
    }

    bool contains(const KaucherInterval& other) const
    {
        if (isEmpty() || other.isEmpty()) return false;
        if (isProper())
        {
            return other.start >= start && other.end <= end;
        }
        return other.end >= start && other.start <= end;
    }

    bool intersects(const KaucherInterval& other) const
    {
        if (isEmpty() || other.isEmpty()) return false;
        if (isProper() && other.isProper())
        {
            return !(end < other.start || start > other.end);
        }
        KaucherInterval p1 = proper();
        KaucherInterval p2 = other.proper();
        return !(p1.end < p2.start || p1.start > p2.end);
    }

    KaucherInterval operator+(const KaucherInterval& other) const
    {
        if (isEmpty() || other.isEmpty()) return KaucherInterval();
        return KaucherInterval(start + other.start, end + other.end);
    }

    KaucherInterval operator-(const KaucherInterval& other) const
    {
        if (isEmpty() || other.isEmpty()) return KaucherInterval();
        return KaucherInterval(start - other.end, end - other.start);
    }

    KaucherInterval algebraicDiff(const KaucherInterval& other) const
    {
        if (isEmpty() || other.isEmpty()) return KaucherInterval();
        return KaucherInterval(start - other.start, end - other.end);
    }

    KaucherInterval operator*(const KaucherInterval& other) const
    {
        if (isEmpty() || other.isEmpty()) return KaucherInterval();

        double a = start, b = end;
        double c = other.start, d = other.end;

        double v1 = a * c;
        double v2 = a * d;
        double v3 = b * c;
        double v4 = b * d;

        double new_start = std::min({v1, v2, v3, v4});
        double new_end = std::max({v1, v2, v3, v4});

        return KaucherInterval(new_start, new_end);
    }

    KaucherInterval operator/(const KaucherInterval& other) const
    {
        if (isEmpty() || other.isEmpty())
            return KaucherInterval();

        if (other.contains(0.0))
        {
            throw std::runtime_error("Interval division by interval containing zero");
        }

        double a = start, b = end;
        double c = other.start, d = other.end;

        KaucherInterval inv;
        if (c > 0)
        {
            inv = KaucherInterval(1.0/d, 1.0/c);
        }
        else if (d < 0)
        {
            inv = KaucherInterval(1.0/c, 1.0/d);
        }
        else
        {
            double c_abs = std::abs(c);
            double d_abs = std::abs(d);
            double max_val = std::max(c_abs, d_abs);
            double min_val = std::min(c_abs, d_abs);
            inv = KaucherInterval(-max_val, -min_val);
        }

        return (*this) * inv;
    }

    KaucherInterval operator*(double s) const
    {
        if (isEmpty()) return KaucherInterval();
        if (s >= 0)
            return KaucherInterval(start * s, end * s);
        return KaucherInterval(end * s, start * s);
    }

    KaucherInterval operator/(double s) const
    {
        if (s == 0)
            throw std::runtime_error("Division by zero");
        return (*this) * (1.0 / s);
    }

    KaucherInterval operator-() const
    {
        if (isEmpty()) return KaucherInterval();
        return KaucherInterval(-end, -start);
    }

    KaucherInterval& operator+=(const KaucherInterval& other)
    {
        if (!isEmpty() && !other.isEmpty())
        {
            start += other.start;
            end += other.end;
        }
        return *this;
    }

    KaucherInterval& operator-=(const KaucherInterval& other)
    {
        if (!isEmpty() && !other.isEmpty())
        {
            start -= other.end;
            end -= other.start;
        }
        return *this;
    }

    KaucherInterval& operator*=(const KaucherInterval& other)
    {
        *this = (*this) * other;
        return *this;
    }

    KaucherInterval& operator*=(double s)
    {
        if (!isEmpty())
        {
            if (s >= 0)
            {
                start *= s;
                end *= s;
            }
            else
            {
                std::swap(start, end);
                start *= s;
                end *= s;
            }
        }
        return *this;
    }

    bool isEqual(const KaucherInterval& other, const double tol) const
    {
        return std::abs(start - other.start) < tol && std::abs(end - other.end) < tol;
    }

    static KaucherInterval intersection(const KaucherInterval& a, const KaucherInterval& b)
    {
        return a.meet(b);
    }

    static KaucherInterval hull(const KaucherInterval& a, const KaucherInterval& b)
    {
        return a.join(b);
    }

    static KaucherInterval pow(const KaucherInterval& x, int n)
    {
        if (x.isEmpty()) return KaucherInterval();

        if (n == 0)
            return KaucherInterval(1.0, 1.0);

        if (n < 0)
            return KaucherInterval::pow(x.proper(), -n).inverse();

        if (n == 1)
            return x;

        if (n == 2)
        {
            double a = x.start, b = x.end;
            if (a >= 0)
                return KaucherInterval(a * a, b * b);
            if (b <= 0)
                return KaucherInterval(b * b, a * a);
            return KaucherInterval(0.0, std::max(a * a, b * b));
        }

        KaucherInterval result = x;
        for (int i = 1; i < n; ++i)
            result = result * x;
        return result;
    }

    KaucherInterval inverse() const
    {
        if (isEmpty())
            return KaucherInterval();

        if (contains(0.0))
            throw std::runtime_error("Cannot compute inverse of interval containing zero");

        double a = start, b = end;
        if (a > 0)
            return KaucherInterval(1.0/b, 1.0/a);
        if (b < 0)
            return KaucherInterval(1.0/b, 1.0/a);
        double a_abs = std::abs(a);
        double b_abs = std::abs(b);
        double max_val = std::max(a_abs, b_abs);
        double min_val = std::min(a_abs, b_abs);
        return KaucherInterval(-max_val, -min_val);
    }
};

inline KaucherInterval operator+(double s, const KaucherInterval& x)
{
    return KaucherInterval(s) + x;
}

inline KaucherInterval operator-(double s, const KaucherInterval& x)
{
    return KaucherInterval(s) - x;
}

inline KaucherInterval operator*(double s, const KaucherInterval& x)
{
    return x * s;
}

inline std::ostream& operator<<(std::ostream& os, const KaucherInterval& di)
{
    os << "[" << di.start << ", " << di.end << "]";
    if (di.isImproper())
        os << " (Improper)";
    return os;
}

namespace interval
{
    inline KaucherInterval abs(const KaucherInterval& x)
    {
        if (x.isEmpty()) return KaucherInterval();
        if (x.isProper())
        {
            if (x.start >= 0) return x;
            if (x.end <= 0) return KaucherInterval(-x.end, -x.start);
            return KaucherInterval(0.0, std::max(-x.start, x.end));
        }
        double a = std::abs(x.start);
        double b = std::abs(x.end);
        return KaucherInterval(std::min(a, b), std::max(a, b));
    }

    inline KaucherInterval sqr(const KaucherInterval& x)
    {
        return x * x;
    }

    inline KaucherInterval sqrt(const KaucherInterval& x)
    {
        if (x.isEmpty()) return KaucherInterval();
        if (x.isProper())
        {
            double a = x.start, b = x.end;
            if (b < 0)
                throw std::runtime_error("Square root of negative interval");
            double new_start = (a > 0) ? std::sqrt(a) : 0.0;
            double new_end = std::sqrt(b);
            return KaucherInterval(new_start, new_end);
        }
        double a = std::sqrt(std::abs(x.start));
        double b = std::sqrt(std::abs(x.end));
        return KaucherInterval(std::min(a, b), std::max(a, b));
    }

    inline KaucherInterval exp(const KaucherInterval& x)
    {
        if (x.isEmpty()) return KaucherInterval();
        return KaucherInterval(std::exp(x.lower()), std::exp(x.upper()));
    }

    inline KaucherInterval log(const KaucherInterval& x)
    {
        if (x.isEmpty()) return KaucherInterval();
        if (x.upper() <= 0)
            throw std::runtime_error("Log of non-positive interval");
        double a = (x.lower() > 0) ? std::log(x.lower()) : -std::numeric_limits<double>::infinity();
        double b = std::log(x.upper());
        return KaucherInterval(a, b);
    }

    inline KaucherInterval log10(const KaucherInterval& x)
    {
        if (x.isEmpty()) return KaucherInterval();
        if (x.upper() <= 0)
            throw std::runtime_error("Log10 of non-positive interval");
        double a = (x.lower() > 0) ? std::log10(x.lower()) : -std::numeric_limits<double>::infinity();
        double b = std::log10(x.upper());
        return KaucherInterval(a, b);
    }

    inline KaucherInterval sin(const KaucherInterval& x)
    {
        if (x.isEmpty()) return KaucherInterval();

        double a = x.lower();
        double b = x.upper();
        double period = 2.0 * M_PI;

        a = std::fmod(a, period);
        b = std::fmod(b, period);

        while (a < 0) { a += period; b += period; }

        double a_mod = a;
        double b_mod = b;

        double pi_2 = M_PI / 2.0;
        double pi_3 = 3.0 * M_PI / 2.0;

        if (a_mod <= pi_2 && b_mod >= pi_2)
            a_mod = -1.0;
        if (a_mod <= pi_3 && b_mod >= pi_3)
            b_mod = 1.0;

        double min_val = -1.0;
        double max_val = 1.0;

        if (b_mod - a_mod >= period)
        {
            return KaucherInterval(-1.0, 1.0);
        }

        auto eval_sin = [](double t) { return std::sin(t); };

        double v_a = eval_sin(a_mod);
        double v_b = eval_sin(b_mod);

        min_val = std::min({v_a, v_b, min_val});
        max_val = std::max({v_a, v_b, max_val});

        return KaucherInterval(min_val, max_val);
    }

    inline KaucherInterval cos(const KaucherInterval& x)
    {
        if (x.isEmpty()) return KaucherInterval();

        double a = x.lower();
        double b = x.upper();
        double period = 2.0 * M_PI;

        a = std::fmod(a, period);
        b = std::fmod(b, period);

        while (a < 0) { a += period; b += period; }

        double pi = M_PI;

        if (a <= pi && b >= pi)
            b = 1.0;

        double min_val = -1.0;
        double max_val = 1.0;

        if (b - a >= period)
        {
            return KaucherInterval(-1.0, 1.0);
        }

        auto eval_cos = [](double t) { return std::cos(t); };

        double v_a = eval_cos(a);
        double v_b = eval_cos(b);

        min_val = std::min({v_a, v_b, min_val});
        max_val = std::max({v_a, v_b, max_val});

        return KaucherInterval(min_val, max_val);
    }

    inline KaucherInterval tan(const KaucherInterval& x)
    {
        if (x.isEmpty()) return KaucherInterval();

        double a = x.lower();
        double b = x.upper();
        double pi = M_PI;

        double a_mod = std::fmod(a, pi);
        double b_mod = std::fmod(b, pi);

        while (a_mod < -pi/2) { a_mod += pi; b_mod += pi; }
        while (a_mod > pi/2) { a_mod -= pi; b_mod -= pi; }

        if (a_mod < pi/2 && b_mod > pi/2)
        {
            return KaucherInterval(-std::numeric_limits<double>::infinity(),
                                         std::numeric_limits<double>::infinity());
        }

        double v_a = std::tan(a_mod);
        double v_b = std::tan(b_mod);

        return KaucherInterval(std::min(v_a, v_b), std::max(v_a, v_b));
    }

    inline KaucherInterval asin(const KaucherInterval& x)
    {
        if (x.isEmpty()) return KaucherInterval();
        if (x.lower() < -1.0 || x.upper() > 1.0)
            throw std::runtime_error("Asin argument out of [-1, 1]");
        return KaucherInterval(std::asin(std::max(-1.0, x.lower())),
                                     std::asin(std::min(1.0, x.upper())));
    }

    inline KaucherInterval acos(const KaucherInterval& x)
    {
        if (x.isEmpty()) return KaucherInterval();
        if (x.lower() < -1.0 || x.upper() > 1.0)
            throw std::runtime_error("Acos argument out of [-1, 1]");
        return KaucherInterval(std::acos(std::min(1.0, x.upper())),
                                     std::acos(std::max(-1.0, x.lower())));
    }

    inline KaucherInterval atan(const KaucherInterval& x)
    {
        if (x.isEmpty()) return KaucherInterval();
        return KaucherInterval(std::atan(x.lower()), std::atan(x.upper()));
    }

    inline KaucherInterval sinh(const KaucherInterval& x)
    {
        if (x.isEmpty()) return KaucherInterval();
        return KaucherInterval(std::sinh(x.lower()), std::sinh(x.upper()));
    }

    inline KaucherInterval cosh(const KaucherInterval& x)
    {
        if (x.isEmpty()) return KaucherInterval();
        double a = std::cosh(x.lower());
        double b = std::cosh(x.upper());
        double min_val = std::min(a, b);
        if (x.lower() <= 0 && x.upper() >= 0)
            min_val = 1.0;
        return KaucherInterval(min_val, std::max(a, b));
    }

    inline KaucherInterval tanh(const KaucherInterval& x)
    {
        if (x.isEmpty()) return KaucherInterval();
        return KaucherInterval(std::tanh(x.lower()), std::tanh(x.upper()));
    }

    inline KaucherInterval pow(const KaucherInterval& x, double p)
    {
        if (x.isEmpty()) return KaucherInterval();

        if (p == 0.0)
            return KaucherInterval(1.0, 1.0);

        if (p == 1.0)
            return x;

        if (p == 2.0)
            return x * x;

        if (p == 0.5)
            return interval::sqrt(x);

        if (x.isProper() && x.start >= 0)
        {
            return KaucherInterval(std::pow(x.start, p), std::pow(x.end, p));
        }

        if (x.isProper() && x.end <= 0 && p == std::floor(p))
        {
            double a = std::pow(-x.end, p);
            double b = std::pow(-x.start, p);
            if (std::fmod(p, 2.0) == 1.0)
                return KaucherInterval(-b, -a);
            return KaucherInterval(a, b);
        }

        double a = x.lower();
        double b = x.upper();
        double v1 = std::pow(a, p);
        double v2 = std::pow(b, p);
        double min_val = std::min(v1, v2);
        double max_val = std::max(v1, v2);

        if (a < 0 && b > 0)
        {
            double v3 = std::pow(-a, p);
            min_val = std::min({min_val, v3});
            max_val = std::max({max_val, v3});
        }

        return KaucherInterval(min_val, max_val);
    }
}
} // namespace ik
