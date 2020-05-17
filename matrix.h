#include <vector>
#include <stdexcept>
#include <iostream>
#include <cmath>

using std::vector;

template <typename T>
vector<T>& operator+=(vector<T>& v1, const vector<T>& v2) {
    if (v1.size() != v2.size())
        throw std::out_of_range{"operator+=: sizes do not match"};

    for (size_t i = 0; i < v1.size(); ++i)
        v1[i] += v2[i];
    return v1;
}

template <typename T>
vector<T> operator+(vector<T> v1, const vector<T>& v2) {
    v1 += v2;
    return v1;
}

template <typename T>
vector<T>& operator-=(vector<T>& v1, const vector<T>& v2) {
    if (v1.size() != v2.size())
        throw std::out_of_range{"operator-=: sizes do not match"};

    for (size_t i = 0; i < v1.size(); ++i)
        v1[i] -= v2[i];
    return v1;
}

template <typename T>
vector<T> operator-(vector<T> v1, const vector<T>& v2) {
    v1 -= v2;
    return v1;
}

template <typename T>
vector<T>& operator*=(vector<T>& v, const T scalar) {
    for (size_t i = 0; i < v.size(); ++i)
        v[i] *= scalar;
    return v;
}

template <typename T>
vector<T>& operator*=(const T scalar, vector<T>& v) {
    return v * scalar;
}

template <typename T>
vector<T> operator*(const T scalar, vector<T> v) {
    v *= scalar;
    return v;
}

// template <typename T>
// vector<vector<T>>& operator*=(vector<vector<T>>& v, vector<vector<T>> v2) {
// }

template <typename T>
vector<vector<T>> operator*(const vector<vector<T>>& lhs, const vector<vector<T>>& rhs) {
    if (lhs.size() && lhs[0].size() != rhs.size())
        throw "operator*: matrices' sizes do not match";
    if (lhs.empty() || lhs[0].empty() || rhs.empty() || rhs[0].empty())
        throw "operator*: degenerate matrix";
    vector<vector<T>> ans(lhs.size(), vector<T>(rhs[0].size(), 0));
    for (size_t i = 0; i < lhs.size(); ++i) {
        for (size_t j = 0; j < rhs[0].size(); ++j) {
            for (size_t k = 0; k < rhs.size(); ++k) {
                ans[i][j] += lhs[i][k] * rhs[k][j];
            }
        }
    }
    return ans;
}

template <typename T>
double euclideanNorm(const vector<T>& v) {
    double ans = 0;
    for (size_t i = 0; i < v.size(); ++i) {
        ans += v[i] * v[i];
    }
    ans /= v.size(); // yes, I know that there isn't such thing in euclidean norm
    return sqrt(ans);
}

template <typename T>
double manhettenNorm(const vector<T>& v) {
    double ans = 0;
    for (size_t i = 0; i < v.size(); ++i) {
        ans += abs(v[i]);
    }
    return ans;
}

template <typename T>
vector<T> operator-(const vector<T>& v) {
    return -1.0 * v;
}

template <typename T>
std::ostream& operator<<(std::ostream& os, const vector<T>& v) {
    os << "{";
    for (size_t i = 0; i < v.size(); ++i) {
        os << v[i];
        if (i + 1 != v.size())
            os << ", ";
    }
    os << "}";
    return os;
}

template <typename T>
std::ostream& operator<<(std::ostream& os, const vector<vector<T>>& v) {
    for (size_t i = 0; i < v.size(); ++i) {
        for (size_t j = 0; j < v[i].size(); ++j) {
            os << v[i][j] << "\t";
        }
        os << "\n";
    }
    return os;
}