#include <iostream>
#include <cmath>
#include <utility>
#include <algorithm>
#include "matrix.h"

class IRK {
private:
    const int stages_;
    const vector<vector<double>> A_;
    const vector<double> b_;
    const vector<double> c_;

    double Rtol_;
    double Atol_;
    int dim_;
    vector<double> (*f_)(double, vector<double>);
    vector<vector<double>> (*Jacobian_)(vector<double>);

    double tn_;
    vector<double> yn_;
    double h_;

    const double EPS_ = 1e-9;
    vector<vector<double>> LU_;
    vector<double> P_;
    double k_;

public:
    IRK(int stages, const vector<vector<double>> A, const vector<double> b,
        const vector<double> c, double Rtol, double Atol, int dim,
        vector<double> (*f)(double, vector<double>),
        vector<vector<double>> (*Jacobian)(vector<double>))
            : stages_(stages)
            , A_(A)
            , b_(b)
            , c_(c)
            , Rtol_(Rtol)
            , Atol_(Atol)
            , dim_(dim)
            , f_(f)
            , Jacobian_(Jacobian) {
        LU_.resize(dim_ * stages_, vector<double>(dim_ * stages_));
        P_.resize(dim_ * stages_);
        k_ = 0.05;
    }

    vector<double> LUPSolve(vector<double> rhs) {
        vector<double> x(rhs.size());
        for (int i = 0; i < rhs.size(); ++i) {
            x[i] = rhs[P_[i]];
            for (int j = 0; j < i; ++j) {
                x[i] -= LU_[i][j] * x[j];
            }
        }

        for (int i = static_cast<int>(rhs.size()) - 1; i >= 0; --i) {
            for (int j = i + 1; j < rhs.size(); ++j) {
                x[i] -= LU_[i][j] * x[j];
            }
            x[i] /= LU_[i][i];
            if (std::isnan(x[i])) {
                throw "nan\n";
            }
        }

        return x;
    }

    void LUPDecompose() {
        for (int i = 0; i < LU_.size(); ++i) {
            P_[i] = i;
        }

        for (int col = 0; col < LU_.size(); ++col) {
            double max_el = 0.0;
            double abs_el;
            int col_max = col;

            for (int k = col; k < LU_.size(); ++k) {
                if ((abs_el = fabs(LU_[k][col])) > max_el) { 
                    max_el = abs_el;
                    col_max = k;
                }
            }

            if (max_el < EPS_) {
                throw "LUPDecompose: matrix is degenerate";
            }

            if (col_max != col) {
                std::swap(P_[col], P_[col_max]);
                std::swap(LU_[col], LU_[col_max]);
            }

            for (int row = col + 1; row < LU_.size(); ++row) {
                LU_[row][col] /= LU_[col][col];
                for (int k = col + 1; k < LU_.size(); ++k) {
                    LU_[row][k] -= LU_[row][col] * LU_[col][k];
                }
            }
        }
    }

    vector<double> NewtonStep(vector<double> zk) {
        vector<vector<double>> zs(stages_);
        for (int i = 0; i < stages_; ++i) {
            zs[i] = vector<double>(zk.begin() + i * dim_, zk.begin() + (i + 1) * dim_);
        }
        vector<double> z;
        for (int i = 0; i < stages_; ++i) {
            vector<double> zi(dim_, 0);
            for (int l = 0; l < stages_; ++l) {
                zi += h_ * A_[i][l] * f_(tn_ + c_[l] * h_, yn_ + zs[l]);
            }
            z.insert(z.end(), zi.begin(), zi.end());
        }
        return LUPSolve(-zk + z);
    }

    void recalcLU(vector<double> point) {
        auto J = Jacobian_(point);
        for (int i = 0; i < stages_; ++i) {
            for (int j = 0; j < stages_; ++j) {
                for (int k = 0; k < dim_; ++k) {
                    for (int l = 0; l < dim_; ++l)
                        LU_[i * dim_ + k][j * dim_ + l] = -h_ * A_[i][j] * J[k][l];
                    if (i == j)
                        LU_[i * dim_ + k][j * dim_ + k] += 1;
                }
            }
        }
        LUPDecompose();
    }

    void properLU(vector<double> zk) {
        vector<vector<double>> zs(stages_);
        for (int i = 0; i < stages_; ++i) {
            zs[i] = vector<double>(zk.begin() + i * dim_, zk.begin() + (i + 1) * dim_);
        }
        for (int column = 0; column < stages_; ++column) {
            auto J = Jacobian_(yn_ + zs[column]);
            for (int row = 0; row < stages_; ++row) {
                for (int k = 0; k < dim_; ++k) {
                    for (int l = 0; l < dim_; ++l) {
                        LU_[row * dim_ + k][column * dim_ + l] = 
                            -h_ * A_[row][column] * J[k][l];
                    }
                    if (row == column) {
                        LU_[row * dim_ + k][column * dim_ + k] += 1;
                    }
                }
            }
        }
        LUPDecompose();
    }

    // Solving ordinary differential equations II 
    // IV.8 Simplified Newton Iterations
    int RKStep() {
        vector<double> zk(dim_ * stages_, 0);
        recalcLU(yn_);
        int newtonIter = 0;
        const int newtonIterMax = 10;
        double prevNorm, norm, convRate = 1.0;
        while (newtonIter <= 1 || (newtonIter < newtonIterMax && convRate < 1.0 &&
                pow(convRate, newtonIterMax - newtonIter) / (1.0 - convRate) * norm <= k_ * Rtol_ && 
                convRate / (1.0 - convRate) * norm > k_ * Rtol_)) {
            vector<double> delta = NewtonStep(zk);
            norm = euclideanNorm(delta);
            convRate = norm / prevNorm;
            prevNorm = norm;
            zk += delta;
            ++newtonIter;
        }
        if (convRate >= 1.0 ||
                pow(convRate, newtonIterMax - newtonIter) / (1.0 - convRate) * norm > k_ * Rtol_) {
            return 1;
        }
        for (int i = 0; i < stages_; ++i) {
            vector<double> tmp(zk.begin() + i * dim_, zk.begin() + (i + 1) * dim_);
            yn_ += h_ * b_[i] * f_(tn_ + c_[i] * h_, yn_ + tmp);
        }
        tn_ += h_;
        return 0;
    }

    double error(const vector<double>& prev_yn, const vector<double>& approx) {
        double ans = 0;
        for (size_t i = 0; i < prev_yn.size(); ++i) {
            ans += pow((yn_[i] - approx[i]) / (Atol_ + std::max(prev_yn[i], yn_[i]) * Rtol_), 2);
        }
        ans /= prev_yn.size();
        return sqrt(ans);
    }

    vector<double> solve(double from, double to, vector<double> y0, double firstStep) {
        int rejectedRk = 0;
        int rejectedNewton = 0;
        int rkIters = 0;
        const double fac = 0.7;
        const double facmin = 1e-6;
        double facmax = 2;
        tn_ = from;
        yn_ = y0;
        h_ = std::min(facmax, std::max(facmin, firstStep));
        while (tn_ < to) {
            ++rkIters;

            h_ = std::max(h_, facmin);
            h_ = std::min(h_, to - tn_);
            
            vector<double> prev_yn = yn_;
            double prev_tn = tn_;

            vector<double> small_step = yn_;
            if (RKStep() || RKStep()) {
                yn_ = std::move(prev_yn);
                tn_ = prev_tn;
                if (h_ <= facmin + EPS_) {
                    throw "RK does not converge";
                }
                h_ /= 2;
                ++rejectedNewton;
                continue;
            }
            std::swap(small_step, yn_);
            tn_ = prev_tn;
            
            h_ *= 2;
            if (RKStep()) {
                yn_ = std::move(prev_yn);
                tn_ = prev_tn;
                h_ /= 2;
                if (h_ <= facmin + EPS_) {
                    throw "RK does not converge";
                }
                h_ /= 2;
                ++rejectedNewton;
                continue;
            }
            h_ /= 2;
            
            double err = error(prev_yn, small_step);
            double facopt = fac * pow(1 / err, 1.0 / (stages_ + 1));
            h_ *= std::min(facmax, std::max(facmin, facopt));
            facmax = 2.0;
            if (err > 1.0) {
                yn_ = std::move(prev_yn);
                tn_ = prev_tn;
                ++rejectedRk;
                // it is advisable to put facmax = 1 in the steps right after step-rejection
                facmax = 1.0;
            }
        }
        // std::cout << "rejected RK method steps: " << rejectedRk << "\n";
        // std::cout << "rejected Newton method steps: " << rejectedNewton << "\n";
        // std::cout << "RK steps: " << rkIters << "\n";
        return yn_;
    }
};

class ImplicitMidpointRule : public IRK {
public:
    ImplicitMidpointRule(double Rtol, double Atol, int dim, 
        vector<double> (*f)(double, vector<double>), 
        vector<vector<double>> (*Jacobian)(vector<double>)) 
        : IRK(1, {{0.5}}, {1.0}, {0.5}, Rtol, Atol, dim, f, Jacobian) {}
};

class DIRK3 : public IRK {
public:
    DIRK3(double Rtol, double Atol, int dim,
        vector<double> (*f)(double, vector<double>), 
        vector<vector<double>> (*Jacobian)(vector<double>)) 
        : IRK(2, 
              {
                {0.5 + 1 / (2 * sqrt(3.0)), 0},
                {-1 / sqrt(3.0), 0.5 + 1 / (2 * sqrt(3.0))}
              },
              {0.5, 0.5}, 
              {0.5 + 1 / (2 * sqrt(3.0)), 0.5 - 1 / (2 * sqrt(3.0))}, 
              Rtol, Atol, dim, f, Jacobian) {
        }
};

// "RADAU5, one of the most popular codes for solving DAEs, 
// is based on the 3-stage Radau method."
class Radau5 : public IRK {
public:
    Radau5(double Rtol, double Atol, int dim, 
        vector<double> (*f)(double, vector<double>), 
        vector<vector<double>> (*Jacobian)(vector<double>)) 
        : IRK(3,
              {
                {1.0 / 9.0, (-1.0 - sqrt(6.0)) / 18.0, (-1.0 + sqrt(6.0)) / 18.0},
                {1.0 / 9.0, (88.0 + 7.0 * sqrt(6.0)) / 360.0, (88.0 - 43.0 * sqrt(6.0)) / 360.0},
                {1.0 / 9.0, (88.0 + 43.0 * sqrt(6.0)) / 360.0, (88.0 - 7.0 * sqrt(6.0)) / 360.0}
              }, 
              {1.0 / 9.0, (16.0 + sqrt(6.0)) / 36.0, (16.0 - sqrt(6.0)) / 36.0},
              {0.0, (6.0 - sqrt(6.0)) / 10.0, (6.0 + sqrt(6.0)) / 10.0},
              Rtol, Atol, dim, f, Jacobian) {}
};
