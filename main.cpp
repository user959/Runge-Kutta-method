#include <vector>
#include <iostream>
#include "runge_kutta.h"

using namespace std;

namespace B1 {
    // B1 has complex eigenvalues, and some of its components
    // decay rapidly while others decay slowly
    const int dim = 4;

    vector<double> solution(double t) {
        return {exp(-t) * cos(10 * t), 
                -10 * exp(-t) * sin(10 * t),
                exp(-100 * t) * cos(100 * t),
                -100 * exp(-100 * t) * sin(100 * t)};
    }

    vector<double> f(double t, vector<double> point) {
        vector<vector<double>> coef = {
            {-1, 1, 0, 0},
            {-100, -1, 0, 0},
            {0, 0, -100, 1},
            {0, 0, -10000, -100}
        };
        vector<double> ans(dim, 0);
        for (int i = 0; i < dim; ++i) {
            for (int j = 0; j < dim; ++j) {
                ans[i] += coef[i][j] * point[j];
            }
        }
        return ans;
    }

    vector<vector<double>> Jacobian(vector<double> point) {
        vector<vector<double>> ans = {
            {-1, 1, 0, 0}, 
            {-100, -1, 0, 0},
            {0, 0, -100, 1},
            {0, 0, -10000, -100}
        };
        return ans;
    }
};

namespace B5 {
    // linear
    const int dim = 6;

    vector<double> solution(double t) {
        return {exp(-10 * t) * sin(100 * t) + exp(-10 * t) * cos(100 * t),
                -exp(-10 * t) * sin(100 * t) + exp(-10 * t) * cos(100 * t),
                exp(-4 * t),
                exp(-t),
                exp(-0.5 * t),
                exp(-0.1 * t)};
    }

    vector<double> f(double t, vector<double> point) {
        const double coef[dim][dim] = {
            {-10, 100, 0, 0, 0, 0},
            {-100, -10, 0, 0, 0, 0},
            {0, 0, -4, 0, 0, 0},
            {0, 0, 0, -1, 0, 0},
            {0, 0, 0, 0, -0.5, 0},
            {0, 0, 0, 0, 0, -0.1}
        };
        vector<double> ans(dim, 0);
        for (int i = 0; i < dim; ++i) {
            for (int j = 0; j < dim; ++j) {
                ans[i] += coef[i][j] * point[j];
            }
        }
        return ans;
    }

    vector<vector<double>> Jacobian(vector<double> p) {
        vector<vector<double>> ans = {
            {-10, 100, 0, 0, 0, 0},
            {-100, -10, 0, 0, 0, 0},
            {0, 0, -4, 0, 0, 0},
            {0, 0, 0, -1, 0, 0},
            {0, 0, 0, 0, -0.5, 0},
            {0, 0, 0, 0, 0, -0.1}
        };
        return ans;
    }
};

namespace C1 {
    // shows nonlinear coupling from transient components to smooth component
    const int dim = 4;
    vector<double> y0 = {1, 1, 1, 1};

    vector<double> f(double t, vector<double> p) {
        vector<double> ans = {
            -p[0] + p[1] * p[1] + p[2] * p[2] + p[3] * p[3],
            -10 * p[1] + 10 * (p[2] * p[2] + p[3] * p[3]),
            -40 * p[2] + 40 * p[3] * p[3],
            -100 * p[3] + 2
        };
        return ans;
    }

    vector<vector<double>> Jacobian(vector<double> p) {
        vector<vector<double>> ans = {
            {-1, 2 * p[1], 2 * p[2], 2 * p[3]},
            {0, -10, 20 * p[2], 20 * p[3]},
            {0, 0, -40, 80 * p[3]},
            {0, 0, 0, -100},
        };
        return ans;
    }
};

namespace ProbablyNotStiff {
    const int dim = 2;

    vector<double> solution(double t) {
        return {2 * exp(2 * t) + exp(-t), 
                exp(2 * t) - exp(-t)};
    }

    vector<double> f(double t, vector<double> point) {
        const double coef[dim][dim] = {
            {1, 2},
            {1, 0}
        };
        vector<double> ans(dim, 0);
        for (int i = 0; i < dim; ++i) {
            for (int j = 0; j < dim; ++j) {
                ans[i] += coef[i][j] * point[j];
            }
        }
        return ans;
    }

    vector<vector<double>> Jacobian(vector<double> point) {
        vector<vector<double>> ans(dim, vector<double>(dim));
        ans[0][0] = point[0];
        ans[0][1] = 2 * point[1];
        ans[1][0] = point[0];
        ans[1][1] = 0;
        return ans;
    }
}

int main() {
    using namespace B5;

    double left = 0;
    double right = 10;
    double firstStep = 0.001;

    double Rtol = 0.01;
    double Atol = 0.00001;
    Radau5 radau(Rtol, Atol, dim, f, Jacobian);
    DIRK3 dirk(Rtol, Atol, dim, f, Jacobian);
    ImplicitMidpointRule mid(Rtol, Atol, dim, f, Jacobian);
    vector<double> radau_ans, mid_ans, dirk_ans;
    try {
        radau_ans = radau.solve(left, right, solution(left), firstStep);
    } catch (char const *s) {
        cout << s << "\n";
    }
    try {
        mid_ans = mid.solve(left, right, solution(left), firstStep);
    } catch (char const *s) {
        cout << s << "\n";
    }
    try {
        dirk_ans = dirk.solve(left, right, solution(left), firstStep);
    } catch (char const *s) {
        cout << s << "\n";
    }
    cout << "mid: " << mid_ans - solution(left) << "\n";
    cout << "radau: " << radau_ans - solution(left) << "\n";
    cout << "dirk: " << dirk_ans - solution(left) << "\n";
    cout << "real: " << solution(right) - solution(left) << "\n";
}