#pragma once

#include <array>
#include <cmath>
#include <iostream>
#include <algorithm>
#include <stdexcept>
#include <Eigen/Dense>
#include "KaucherInterval.h"
#include "IntervalVector.h"

namespace ik {

// Fixed-size interval matrix class
template<size_t M, size_t N>
class IntervalMatrix
{
public:
    std::array<std::array<KaucherInterval, N>, M> data;

    static constexpr size_t Rows = M;
    static constexpr size_t Cols = N;

    IntervalMatrix()
    {
        for (size_t i = 0; i < M; ++i)
            for (size_t j = 0; j < N; ++j)
                data[i][j] = KaucherInterval();
    }

    static IntervalMatrix<M, N> zeros()
    {
        IntervalMatrix<M, N> result;
        for (size_t i = 0; i < M; ++i)
            for (size_t j = 0; j < N; ++j)
                result.data[i][j] = KaucherInterval(0.0, 0.0);
        return result;
    }

    static IntervalMatrix<M, N> ones()
    {
        IntervalMatrix<M, N> result;
        for (size_t i = 0; i < M; ++i)
            for (size_t j = 0; j < N; ++j)
                result.data[i][j] = KaucherInterval(1.0, 1.0);
        return result;
    }

    static IntervalMatrix<M, N> identity()
    {
        IntervalMatrix<M, N> result;
        for (size_t i = 0; i < std::min(M, N); ++i)
            result.data[i][i] = KaucherInterval(1.0, 1.0);
        return result;
    }

    static IntervalMatrix<M, N> constants(double c)
    {
        IntervalMatrix<M, N> result;
        for (size_t i = 0; i < M; ++i)
            for (size_t j = 0; j < N; ++j)
                result.data[i][j] = KaucherInterval(c, c);
        return result;
    }

    static IntervalMatrix<M, N> constants(double lower, double upper)
    {
        IntervalMatrix<M, N> result;
        for (size_t i = 0; i < M; ++i)
            for (size_t j = 0; j < N; ++j)
                result.data[i][j] = KaucherInterval(lower, upper);
        return result;
    }

    KaucherInterval& operator()(size_t i, size_t j)
    {
        return data[i][j];
    }

    const KaucherInterval& operator()(size_t i, size_t j) const
    {
        return data[i][j];
    }

    IntervalMatrix<M, N> operator+(const IntervalMatrix<M, N>& other) const
    {
        IntervalMatrix<M, N> result;
        for (size_t i = 0; i < M; ++i)
            for (size_t j = 0; j < N; ++j)
                result.data[i][j] = data[i][j] + other.data[i][j];
        return result;
    }

    IntervalMatrix<M, N> operator-(const IntervalMatrix<M, N>& other) const
    {
        IntervalMatrix<M, N> result;
        for (size_t i = 0; i < M; ++i)
            for (size_t j = 0; j < N; ++j)
                result.data[i][j] = data[i][j] - other.data[i][j];
        return result;
    }

    IntervalMatrix<M, N> operator-() const
    {
        IntervalMatrix<M, N> result;
        for (size_t i = 0; i < M; ++i)
            for (size_t j = 0; j < N; ++j)
                result.data[i][j] = -data[i][j];
        return result;
    }

    IntervalMatrix<M, N> operator*(double s) const
    {
        IntervalMatrix<M, N> result;
        for (size_t i = 0; i < M; ++i)
            for (size_t j = 0; j < N; ++j)
                result.data[i][j] = data[i][j] * s;
        return result;
    }

    friend IntervalMatrix<M, N> operator*(double s, const IntervalMatrix<M, N>& m)
    {
        return m * s;
    }

    IntervalMatrix<M, N> operator/(double s) const
    {
        if (s == 0.0)
            throw std::runtime_error("Division by zero");
        return (*this) * (1.0 / s);
    }

    IntervalMatrix<M, N>& operator+=(const IntervalMatrix<M, N>& other)
    {
        for (size_t i = 0; i < M; ++i)
            for (size_t j = 0; j < N; ++j)
                data[i][j] += other.data[i][j];
        return *this;
    }

    IntervalMatrix<M, N>& operator-=(const IntervalMatrix<M, N>& other)
    {
        for (size_t i = 0; i < M; ++i)
            for (size_t j = 0; j < N; ++j)
                data[i][j] -= other.data[i][j];
        return *this;
    }

    IntervalMatrix<M, N>& operator*=(double s)
    {
        for (size_t i = 0; i < M; ++i)
            for (size_t j = 0; j < N; ++j)
                data[i][j] *= s;
        return *this;
    }

    IntervalMatrix<M, N> intersection(const IntervalMatrix<M, N>& other) const
    {
        IntervalMatrix<M, N> result;
        for (size_t i = 0; i < M; ++i)
            for (size_t j = 0; j < N; ++j)
                result.data[i][j] = data[i][j].meet(other.data[i][j]);
        return result;
    }

    static IntervalMatrix<M, N> hull(const IntervalMatrix<M, N>& a, const IntervalMatrix<M, N>& b)
    {
        IntervalMatrix<M, N> result;
        for (size_t i = 0; i < M; ++i)
            for (size_t j = 0; j < N; ++j)
                result.data[i][j] = a.data[i][j].join(b.data[i][j]);
        return result;
    }

    bool contains(const IntervalMatrix<M, N>& other) const
    {
        for (size_t i = 0; i < M; ++i)
            for (size_t j = 0; j < N; ++j)
                if (!data[i][j].contains(other.data[i][j]))
                    return false;
        return true;
    }

    bool empty() const
    {
        return data[0][0].isEmpty();
    }

    double maxWidth() const
    {
        double max_w = 0.0;
        for (size_t i = 0; i < M; ++i)
            for (size_t j = 0; j < N; ++j)
            {
                double w = data[i][j].magnitude();
                if (w > max_w)
                    max_w = w;
            }
        return max_w;
    }

    IntervalMatrix<M, N> copy() const
    {
        IntervalMatrix<M, N> result;
        for (size_t i = 0; i < M; ++i)
            for (size_t j = 0; j < N; ++j)
                result.data[i][j] = data[i][j];
        return result;
    }

    void print(std::ostream& os = std::cout) const
    {
        os << "[" << M << "x" << N << " Matrix:\n";
        for (size_t i = 0; i < M; ++i)
        {
            os << "  [";
            for (size_t j = 0; j < N; ++j)
            {
                os << data[i][j];
                if (j < N - 1)
                    os << ", ";
            }
            os << "]";
            if (i < M - 1)
                os << "\n";
        }
        os << "]";
    }

    friend std::ostream& operator<<(std::ostream& os, const IntervalMatrix<M, N>& m)
    {
        m.print(os);
        return os;
    }

    Eigen::MatrixXd midpoint() const
    {
        Eigen::MatrixXd result(M, N);
        for (size_t i = 0; i < M; ++i)
            for (size_t j = 0; j < N; ++j)
                result(static_cast<int>(i), static_cast<int>(j)) = data[i][j].middle();
        return result;
    }

    Eigen::MatrixXd lower() const
    {
        Eigen::MatrixXd result(M, N);
        for (size_t i = 0; i < M; ++i)
            for (size_t j = 0; j < N; ++j)
                result(static_cast<int>(i), static_cast<int>(j)) = data[i][j].lower();
        return result;
    }

    Eigen::MatrixXd upper() const
    {
        Eigen::MatrixXd result(M, N);
        for (size_t i = 0; i < M; ++i)
            for (size_t j = 0; j < N; ++j)
                result(static_cast<int>(i), static_cast<int>(j)) = data[i][j].upper();
        return result;
    }

    Eigen::MatrixXd toEigen() const
    {
        return midpoint();
    }

    static IntervalMatrix<M, N> fromEigen(const Eigen::MatrixXd& mat)
    {
        IntervalMatrix<M, N> result;
        for (size_t i = 0; i < M; ++i)
            for (size_t j = 0; j < N; ++j)
                result.data[i][j] = KaucherInterval(mat(static_cast<int>(i), static_cast<int>(j)), mat(static_cast<int>(i), static_cast<int>(j)));
        return result;
    }
};

// Matrix multiplication operator (fixed size)
template<size_t M1, size_t N1, size_t K1>
IntervalMatrix<M1, N1> operator*(const IntervalMatrix<M1, K1>& A, const IntervalMatrix<K1, N1>& B)
{
    IntervalMatrix<M1, N1> result;
    for (size_t i = 0; i < M1; ++i)
        for (size_t j = 0; j < N1; ++j)
        {
            KaucherInterval sum(0.0, 0.0);
            for (size_t k = 0; k < K1; ++k)
                sum = sum + A.data[i][k] * B.data[k][j];
            result.data[i][j] = sum;
        }
    return result;
}

// Matrix-vector multiplication operator (fixed size)
template<size_t M2, size_t N2>
IntervalVector<M2> operator*(const IntervalMatrix<M2, N2>& A, const IntervalVector<N2>& x)
{
    IntervalVector<M2> result;
    for (size_t i = 0; i < M2; ++i)
    {
        KaucherInterval sum(0.0, 0.0);
        for (size_t j = 0; j < N2; ++j)
            sum = sum + A.data[i][j] * x[j];
        result[i] = sum;
    }
    return result;
}

template<size_t M3, size_t N3>
struct IntervalMatrixAbs
{
    static IntervalMatrix<M3, N3> abs(const IntervalMatrix<M3, N3>& m)
    {
        IntervalMatrix<M3, N3> result;
        for (size_t i = 0; i < M3; ++i)
            for (size_t j = 0; j < N3; ++j)
                result.data[i][j] = KaucherInterval(m.data[i][j].absoluteValue());
        return result;
    }
};

template<size_t M4, size_t N4>
IntervalMatrix<M4, N4> matrix_abs(const IntervalMatrix<M4, N4>& m)
{
    return IntervalMatrixAbs<M4, N4>::abs(m);
}

template<size_t N5>
IntervalMatrix<N5, N5> matrix_identity()
{
    return IntervalMatrix<N5, N5>::identity();
}

namespace matrix
{
    template<size_t M, size_t N>
    class FixedMatrixOperations
    {
    public:
        static Eigen::MatrixXd toDoubleMatrix(const IntervalMatrix<M, N>& m)
        {
            return m.midpoint();
        }

        static Eigen::VectorXd vectorToDouble(const IntervalVector<M>& v)
        {
            return v.midpoint();
        }

        static Eigen::MatrixXd multiply(const Eigen::MatrixXd& A, const Eigen::MatrixXd& B)
        {
            return A * B;
        }

        static Eigen::VectorXd multiply(const Eigen::MatrixXd& A, const Eigen::VectorXd& x)
        {
            return A * x;
        }

        static Eigen::MatrixXd inverse(const Eigen::MatrixXd& A)
        {
            return A.inverse();
        }

        static Eigen::MatrixXd transpose(const Eigen::MatrixXd& A)
        {
            return A.transpose();
        }

        static double determinant(const Eigen::MatrixXd& A)
        {
            return A.determinant();
        }

        static bool isInvertible(const Eigen::MatrixXd& A, double tol = 1e-12)
        {
            return std::abs(A.determinant()) > tol;
        }

        static Eigen::MatrixXd identity(size_t n)
        {
            return Eigen::MatrixXd::Identity(static_cast<int>(n), static_cast<int>(n));
        }

        static Eigen::MatrixXd zeros(size_t m, size_t n)
        {
            return Eigen::MatrixXd::Zero(static_cast<int>(m), static_cast<int>(n));
        }
    };
}
} // namespace ik