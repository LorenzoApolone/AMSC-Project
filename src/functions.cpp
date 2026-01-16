/**
 * @file functions_part2.hpp
 * @brief Test functions (benchmark) with Doxygen-style documentation.
 */
#include "interfaces.hpp"

/**
 * @class Sphere
 * @brief Sphere (De Jong) function: f(x) = Σ x_i^2. Convex, separable; global
 * minimum at x* = 0.
 */
class Sphere : public TestFunction {

public:
  /**
   * @brief Construct the Sphere test function.
   * @param dim Problem dimension.
   */
  Sphere(unsigned int dim)
      : TestFunction(dim, "Sphere", std::pair<double, double>{-100.0, 100.0},
                     std::vector<double>(dim, 0.)) {};
  /**
   * @brief Evaluate f(x) = Σ x_i^2.
   * @param x Input vector (size = dim).
   * @return Function value at x.
   */
  double value(const std::vector<double> &x) const override {
    double sum = 0.0;
    for (double xi : x) {
      sum += xi * xi;
    }
    return sum;
  }
};
/**
 * @class Ellipsoid
 * @brief Weighted sphere (ellipsoid): f(x) = Σ (i+1)·x_i^2. Ill-conditioned as
 * i grows.
 */
class Ellipsoid : public TestFunction {

public:
  /**
   * @brief Construct the Ellipsoid test function.
   * @param dim Problem dimension.
   */
  Ellipsoid(unsigned int dim)
      : TestFunction(dim, "Ellipsoid", std::pair<double, double>{-100.0, 100.0},
                     std::vector<double>(dim, 0.)) {};
  /**
   * @brief Evaluate f(x) = Σ (i+1)·x_i^2.
   * @param x Input vector (size = dim).
   * @return Function value at x.
   */
  double value(const std::vector<double> &x) const override {
    double sum = 0.0;
    for (size_t i = 0; i < dim; i++) {
      sum += (i + 1) * x[i] * x[i];
    }
    return sum;
  }
};

/**
 * @class SumOfDiffPowers
 * @brief Sum of different powers: f(x) = Σ |x_i|^{i+2}. Increasing nonlinearity
 * across coordinates.
 */

class SumOfDiffPowers : public TestFunction {

public:
  SumOfDiffPowers(unsigned int dim)
      : TestFunction(dim, "SumOfDiffPowers",
                     std::pair<double, double>{-10.0, 10.0},
                     std::vector<double>(dim, 0.)) {};

  double value(const std::vector<double> &x) const override {
    double sum = 0.0;
    for (size_t i = 0; i < dim; i++) {
      sum += pow(abs(x[i]), i + 2);
    }
    return sum;
  }
};

/**
 * @class QuinticFunction
 * @brief Absolute quintic sum: f(x) = Σ |x_i^5 - 3x_i^4 + 4x_i^3 + 2x_i^2 -
 * 10x_i - 4|.
 */

class QuinticFunction : public TestFunction {

public:
  QuinticFunction(unsigned int dim)
      : TestFunction(dim, "QuinticFunction",
                     std::pair<double, double>{-20.0, 20.0},
                     std::vector<double>(dim, 0.0)) {};

  double value(const std::vector<double> &x) const override {
    double sum = 0.0;
    for (std::size_t i = 0; i < x.size(); ++i) {
      double xi = x[i];
      double term = std::pow(xi, 5) - 3.0 * std::pow(xi, 4) +
                    4.0 * std::pow(xi, 3) + 2.0 * std::pow(xi, 2) - 10.0 * xi -
                    4.0;
      sum += std::fabs(term);
    }
    return sum;
  }
};

/**
 * @class DropWave
 * @brief Drop-Wave: f(x) = 1 - (1 + cos(12·√(Σ x_i^2))) / (2 + 0.5·Σ x_i^2).
 * Multimodal, radially symmetric.
 */

class DropWave : public TestFunction {

public:
  DropWave(unsigned int dim)
      : TestFunction(dim, "DropWave", std::pair<double, double>{-5.12, 5.12},
                     std::vector<double>(dim, 0.)) {};

  double value(const std::vector<double> &x) const override {
    double sum = 1.0;
    double num = 1.;
    double den = 2.;
    double vec_squared = 0.;

    for (size_t i = 0; i < dim; i++) {
      vec_squared += x[i] * x[i];
    }

    num += cos(12 * sqrt(vec_squared));
    den += 0.5 * vec_squared;

    return sum - num / den;
  }
};

/**
 * @class Weierstrass
 * @brief Weierstrass function: highly multimodal, non-differentiable-like
 * behavior via cosine series.
 */
class Weierstrass : public TestFunction {

public:
  Weierstrass(unsigned int dim)
      : TestFunction(dim, "Weierstrass", std::pair<double, double>{-0.5, 0.5},
                     std::vector<double>(dim, 0.)) {};

  double value(const std::vector<double> &x) const override {
    double sum = 0.;
    double a = 0.5;
    unsigned int b = 3;
    unsigned int kmax = 20;

    for (size_t i = 0; i < dim; i++) {
      for (size_t j = 0; j < kmax; j++) {
        sum += pow(a, j) * cos(2 * M_PI * pow(b, j) * (x[i] + 0.5));
      }
    }

    for (size_t i = 0; i < dim; i++) {
      for (size_t j = 0; j < kmax; j++) {
        sum -= dim * pow(a, j) * cos(M_PI * pow(b, j));
      }
    }

    return sum;
  }
};

/**
 * @class Alpine1
 * @brief Alpine #1: f(x) = Σ |x_i·sin(x_i) + 0.1·x_i|. Many local minima.
 */
class Alpine1 : public TestFunction {

public:
  Alpine1(unsigned int dim)
      : TestFunction(dim, "Alpine1", std::pair<double, double>{-10., 10.},
                     std::vector<double>(dim, 0.)) {};

  double value(const std::vector<double> &x) const override {
    double sum = 0.;

    for (size_t i = 0; i < dim; i++) {
      sum += abs(x[i] * sin(x[i]) + 0.1 * x[i]);
    }

    return sum;
  }
};

/**
 * @class Ackley
 * @brief Ackley function: widely used multimodal benchmark with exponential
 * terms.
 */
class Ackley : public TestFunction {

public:
  Ackley(unsigned int dim)
      : TestFunction(dim, "Ackley", std::pair<double, double>{-32.768, 32.768},
                     std::vector<double>(dim, 0.)) {};

  double value(const std::vector<double> &x) const override {
    const double a = 20.0;
    const double b = 0.2;
    const double c = 2.0 * M_PI;
    double sum_sq = 0.0;
    double sum_cos = 0.0;
    for (double xi : x) {
      sum_sq += xi * xi;
      sum_cos += std::cos(c * xi);
    }
    double n = static_cast<double>(x.size());
    double term1 = -a * std::exp(-b * std::sqrt(sum_sq / n));
    double term2 = -std::exp(sum_cos / n);
    return term1 + term2 + a + std::exp(1.0);
  }
};

/**
 * @class Griewank
 * @brief Griewank: f(x) = 1 + Σ(x_i^2/4000) − Π cos(x_i/√(i+1)). Many regularly
 * spaced minima.
 */
class Griewank : public TestFunction {

public:
  Griewank(unsigned int dim)
      : TestFunction(dim, "Griewank", std::pair<double, double>{-100., 100.},
                     std::vector<double>(dim, 0.)) {};

  double value(const std::vector<double> &x) const override {
    double sum = 1.;
    double prod = 1.;

    for (size_t i = 0; i < dim; i++) {
      sum += x[i] * x[i] / 4000.;
      prod *= cos(x[i] / sqrt(i + 1));
    }

    return sum - prod;
  }
};

/**
 * @class Rastrigin
 * @brief Rastrigin: f(x) = 10D + Σ[x_i^2 − 10·cos(2πx_i)]. Highly multimodal,
 * separable.
 */
class Rastrigin : public TestFunction {

public:
  Rastrigin(unsigned int dim)
      : TestFunction(dim, "Rastrigin", std::pair<double, double>{-5.12, 5.12},
                     std::vector<double>(dim, 0.)) {};

  double value(const std::vector<double> &x) const override {
    double sum = 10.0 * dim;
    for (double xi : x) {
      sum += std::pow(xi, 2) - 10.0 * std::cos(2.0 * M_PI * xi);
    }
    return sum;
  }
};

/**
 * @class HappyCat
 * @brief HappyCat (HC): non-separable, with quartic-root ridge term and
 * linear/quadratic coupling.
 */
class HappyCat : public TestFunction {

public:
  HappyCat(unsigned int dim)
      : TestFunction(dim, "HappyCat", std::pair<double, double>{-20., 20.},
                     std::vector<double>(dim, -1.)) {};

  double value(const std::vector<double> &x) const override {
    double sum = -dim;
    double vec_squared = 0.;
    double vec_sum = 0.;

    for (size_t i = 0; i < dim; i++) {
      vec_squared += x[i] * x[i];
      vec_sum += x[i];
    }

    sum += vec_squared;
    sum = pow(abs(sum), 0.25);
    sum += (0.5 * vec_squared + vec_sum) / dim + 0.5;

    return sum;
  }
};

/**
 * @class HGBat
 * @brief HGBat: involves sqrt(|‖x‖_2^4 − (Σ x_i)^2|) + linear/quadratic
 * coupling.
 */
class HGBat : public TestFunction {

public:
  HGBat(unsigned int dim)
      : TestFunction(dim, "HGBat", std::pair<double, double>{-15., 15.},
                     std::vector<double>(dim, -1.)) {};

  double value(const std::vector<double> &x) const override {
    double sum = 0.;
    double vec_squared = 0.;
    double vec_sum = 0.;

    for (size_t i = 0; i < dim; i++) {
      vec_squared += x[i] * x[i];
      vec_sum += x[i];
    }

    sum += sqrt(abs(pow(vec_squared, 2) - pow(vec_sum, 2)));
    sum += (0.5 * vec_squared + vec_sum) / dim + 0.5;

    return sum;
  }
};
/**
 * @class Rosenbrock
 * @brief Rosenbrock “banana” function: narrow curved valley; global min at x*
 * = 1.
 */
class Rosenbrock : public TestFunction {
public:
  Rosenbrock(unsigned int dim)
      : TestFunction(dim, "Rosenbrock", std::pair<double, double>{-10.0, 10.0},
                     std::vector<double>(dim, 1.0)) {}

  double value(const std::vector<double> &x) const override {
    double sum = 0.0;

    for (std::size_t i = 0; i + 1 < dim; ++i) {
      const double xi = x[i];
      const double xip1 = x[i + 1];
      const double t1 = xip1 - xi * xi;
      const double t2 = xi - 1.0;
      sum += 100.0 * t1 * t1 + t2 * t2;
    }
    return sum;
  }
};
/**
 * @class HighCondElliptic
 * @brief High-conditioned elliptic: f(x) = Σ 10^{6·i/(D-1)}·x_i^2, i = 0..D-1.
 */
class HighCondElliptic : public TestFunction {
public:
  HighCondElliptic(unsigned int dim)
      : TestFunction(dim, "HighConditionedElliptic",
                     std::pair<double, double>{-100.0, 100.0},
                     std::vector<double>(dim, 0.0)) {}

  double value(const std::vector<double> &x) const override {
    double sum = 0.0;
    const double exponent_scale = 6.0;
    if (dim == 0)
      return 0.0;
    if (dim == 1)
      return x[0] * x[0];
    for (std::size_t i = 0; i < dim; ++i) {
      double exponent = (static_cast<double>(i) / (dim - 1.0)) * exponent_scale;
      sum += std::pow(10.0, exponent) * x[i] * x[i];
    }
    return sum;
  }
};

/**
 * @class Discus
 * @brief Discus (tablet) function: f(x) = 10^6·x_1^2 + Σ_{i=2..D} x_i^2. Strong
 * axis anisotropy.
 */
class Discus : public TestFunction {
public:
  Discus(unsigned int dim)
      : TestFunction(dim, "Discus", std::pair<double, double>{-100.0, 100.0},
                     std::vector<double>(dim, 0.0)) {}

  double value(const std::vector<double> &x) const override {
    if (x.empty())
      return 0.0;

    double sum = 0.0;
    sum += std::pow(10.0, 6.0) * x[0] * x[0];

    for (std::size_t i = 1; i < dim; ++i) {
      sum += x[i] * x[i];
    }

    return sum;
  }
};

/**
 * @class BentCigar
 * @brief Bent Cigar: f(x) = x_1^2 + 10^6 Σ_{i=2..D} x_i^2. Strong axis
 * anisotropy.
 * @note Global minimum at x* = 0 with f(x*) = 0. Non-separable.
 */
class BentCigar : public TestFunction {
public:
  BentCigar(unsigned int dim)
      : TestFunction(dim, "BentCigar", std::pair<double, double>{-100.0, 100.0},
                     std::vector<double>(dim, 0.0)) {}

  double value(const std::vector<double> &x) const override {
    if (x.empty())
      return 0.0;

    double sum = 0.0;
    sum += x[0] * x[0];
    for (std::size_t i = 1; i < dim; ++i) {
      sum += std::pow(10.0, 6.0) * x[i] * x[i];
    }

    return sum;
  }
};
/**
 * @class PermdbFunc
 * @brief Perm (β) function (decomposable variant).
 *        f(x) = Σ_{i=1..D} [ Σ_{j=1..D} (j^{i} + β) ((x_j/j)^{i} − 1) ]^2
 * @note Global minimum at x*_j = j with f(x*) = 0 (for this parameterization).
 */

class PermdbFunc : public TestFunction {
public:
  PermdbFunc(unsigned int dim)
      : TestFunction(dim, "PermBetaFunction",
                     std::pair<double, double>{-dim, dim}, [dim]() {
                       std::vector<double> v(dim);
                       for (unsigned int i = 0; i < dim; ++i)
                         v[i] = static_cast<double>(i + 1);
                       return v;
                     }()) {}
  double value(const std::vector<double> &x) const override {
    if (x.empty())
      return 0.0;
    double beta = 0.5;

    double sum = 0.0;
    double temp = 0.0;
    for (std::size_t i = 0; i < dim; ++i) {
      temp = 0.0;
      for (std::size_t j = 0; j < dim; ++j) {
        temp += (std::pow(j + 1, i + 1) + beta) *
                (std::pow((x[j] / (j + 1)), i + 1) - 1.0);
      }
      temp *= temp;
      sum += temp;
    }

    return sum;
  }
};
/**
 * @class Schafferf7Func
 * @brief Schaffer F7 (averaged pairwise form).
 *        f(x) = [ (1/(D−1)) Σ_{i=1..D−1} ( √s_i + √s_i·sin^2(50·s_i^{0.2}) )
 * ]^2, with s_i = x_i^2 + x_{i+1}^2.
 * @note Highly multimodal; global minimum at x* = 0.
 */
class Schafferf7Func : public TestFunction {
public:
  Schafferf7Func(unsigned int dim)
      : TestFunction(dim, "SchafferF7",
                     std::pair<double, double>{-100.0, 100.0},
                     std::vector<double>(dim, 0.0)) {}
  double value(const std::vector<double> &x) const override {
    const std::size_t D = x.size();
    if (D < 2)
      return 0.0;

    double acc = 0.0;
    for (std::size_t i = 0; i + 1 < D; ++i) {
      const double xi = x[i];
      const double xip = x[i + 1];

      const double si = std::sqrt(xi * xi + xip * xip);
      const double sqrt_si = std::sqrt(si);
      const double s15 = std::pow(si, 0.2);

      acc += sqrt_si + sqrt_si * std::sin(50.0 * s15) * std::sin(50.0 * s15);
    }

    const double mean = acc * (1.0 / static_cast<double>(D - 1));
    return mean * mean;
  }
};
/**
 * @class ExpSchafferF6
 * @brief Expanded Schaffer F6 with cyclic pairs.
 *        f(x) = Σ_{i=1..D} g(x_i, x_{i+1}), with x_{D+1} = x_1,
 *        g(u,v) = 0.5 + (sin^2(√(u^2+v^2)) − 0.5) / (1 + 0.001(u^2+v^2))^2.
 */
class ExpSchafferF6 : public TestFunction {
public:
  ExpSchafferF6(unsigned int dim)
      : TestFunction(dim, "ExpandedSchafferF6",
                     std::pair<double, double>{-100.0, 100.0},
                     std::vector<double>(dim, 0.0)) {}

  double value(const std::vector<double> &x) const override {
    const std::size_t D = x.size();
    if (D < 2)
      return 0.0;

    double sum = 0.0;
    for (std::size_t i = 0; i + 1 < D; ++i)
      sum += g(x[i], x[i + 1]);

    sum += g(x[D - 1], x[0]);
    return sum;
  }

private:
  static double g(double x, double y) noexcept {
    const double r2 = x * x + y * y;
    const double s = std::sin(std::sqrt(r2));
    const double num = s * s - 0.5;
    const double den = std::pow(1.0 + 0.001 * r2, 2.0);
    return 0.5 + num / den;
  }
};

/**
 * @class RotatedHyper
 * @brief Rotated hyper-ellipsoid (monotone weights): f(x) = Σ_{i=1..D}
 * (D+2−i)·x_i^2.
 * @note Equivalent to a weighted quadratic; global minimum at x* = 0.
 */
class RotatedHyper : public TestFunction {
public:
  RotatedHyper(unsigned int dim)
      : TestFunction(dim, "RotatedHyper-ellipsoidFunction",
                     std::pair<double, double>{-100.0, 100.0},
                     std::vector<double>(dim, 0.0)) {}
  double value(const std::vector<double> &x) const override {
    if (x.empty())
      return 0.0;

    double sum = 0.0;

    for (unsigned int i = 0; i < dim; ++i) {
      sum += (dim + 2 - i) * std::pow(x[i], 2);
    }

    return sum;
  }
};
/**
 * @class Schwefel
 * @brief Schwefel 2.26: f(x) = 418.982887...·D − Σ x_i·sin(√|x_i|).
 * @note Many local minima; global minimum near x_i ≈ 420.968746... for all i.
 */
class Schwefel : public TestFunction {
public:
  Schwefel(unsigned int dim)
      : TestFunction(dim, " Schwefel", std::pair<double, double>{-500.0, 500.0},
                     std::vector<double>(dim, 420.968746359982025)) {}
  double value(const std::vector<double> &x) const override {
    if (x.empty())
      return 0.0;

    double sum = 0.0;
    for (double xi : x)
      sum += xi * sin(sqrt(fabs(xi)));
    return 418.9828872724337 * x.size() - sum;
  }
};
/**
 * @class SumOfDifferentPowers2
 * @brief Sum of different powers with dimension-scaled exponents.
 *        f(x) = √( Σ |x_i|^{ 2 + 4·(i/(D−1)) } ).
 */
class SumOfDifferentPowers2 : public TestFunction {
public:
  SumOfDifferentPowers2(unsigned int dim)
      : TestFunction(dim, "SumOfDifferentPowers",
                     std::pair<double, double>{-10.0, 10.0},
                     std::vector<double>(dim, 0.0)) {}

  double value(const std::vector<double> &x) const override {
    if (x.empty())
      return 0.0;

    double sum = 0.0;
    unsigned int D = x.size();

    for (unsigned int i = 0; i < D; ++i) {
      double exponent =
          2.0 + 4.0 * static_cast<double>(i) / static_cast<double>(D - 1);
      sum += pow(fabs(x[i]), exponent);
    }

    return sqrt(sum);
  }
};
/**
 * @class XinSheYang1
 * @brief Xin-She Yang #1: f(x) = (Σ |x_i|) · exp( − Σ sin(x_i^2) ).
 * @note Multimodal, separable in the exponential argument.
 */
class XinSheYang1 : public TestFunction {
public:
  XinSheYang1(unsigned int dim)
      : TestFunction(dim, "Xin-SheYang's1",
                     std::pair<double, double>{-2.0 * PI(), 2.0 * PI()},
                     std::vector<double>(dim, 0.0)) {}

  double value(const std::vector<double> &x) const override {
    if (x.empty())
      return 0.0;

    double sum_abs = 0.0;
    double sum_sin_sq = 0.0;

    for (double xi : x) {
      sum_abs += std::fabs(xi);
      sum_sin_sq += std::sin(xi * xi);
    }

    return sum_abs * std::exp(-sum_sin_sq);
  }

private:
  static constexpr double PI() { return 3.14159265358979323846; }
};
/**
 * @class Schwefel221
 * @brief Schwefel 2.21: f(x) = max_i |x_i|.
 * @note Non-differentiable at the optimum; global minimum at x* = 0.
 */
class Schwefel221 : public TestFunction {
public:
  Schwefel221(unsigned int dim)
      : TestFunction(dim, "Schwefel2_21",
                     std::pair<double, double>{-100.0, 100.0},
                     std::vector<double>(dim, 0.0)) {}

  double value(const std::vector<double> &x) const override {
    if (x.empty())
      return 0.0;

    double max_abs = 0.0;
    for (double xi : x) {
      double a = std::fabs(xi);
      if (a > max_abs)
        max_abs = a;
    }
    return max_abs;
  }
};
/**
 * @class Schwefel222
 * @brief Schwefel 2.22: f(x) = Σ |x_i| + Π |x_i|.
 * @note Non-separable due to the product term; minimum at x* = 0.
 */
class Schwefel222 : public TestFunction {
public:
  Schwefel222(unsigned int dim)
      : TestFunction(dim, "Schwefel2_22",
                     std::pair<double, double>{-100.0, 100.0},
                     std::vector<double>(dim, 0.0)) {}

  double value(const std::vector<double> &x) const override {
    if (x.empty())
      return 0.0;

    double sum_abs = 0.0;
    double prod_abs = 1.0;

    for (double xi : x) {
      double a = std::fabs(xi);
      sum_abs += a;
      prod_abs *= a;
    }

    return sum_abs + prod_abs;
  }
};
/**
 * @class Salomon
 * @brief Salomon: f(x) = 1 − cos(2π‖x‖_2) + 0.1‖x‖_2.
 * @note Radially symmetric; many concentric local minima.
 */
class Salomon : public TestFunction {
public:
  Salomon(unsigned int dim)
      : TestFunction(dim, "Salomon", std::pair<double, double>{-20.0, 20.0},
                     std::vector<double>(dim, 0.0)) {}

  double value(const std::vector<double> &x) const override {
    if (x.empty())
      return 0.0;

    double sum_sq = 0.0;
    for (double xi : x)
      sum_sq += xi * xi;

    double sqrt_sum = std::sqrt(sum_sq);
    return 1.0 - std::cos(2.0 * PI() * sqrt_sum) + 0.1 * sqrt_sum;
  }

private:
  static constexpr double PI() { return 3.14159265358979323846; }
};
/**
 * @class ModifiedRidge
 * @brief Modified ridge: f(x) = |x_1| + 2·( Σ_{i=2..D} x_i^2 )^{0.1}.
 * @note Non-differentiable at x_1 = 0; ridge-like landscape.
 */
class ModifiedRidge : public TestFunction {
public:
  ModifiedRidge(unsigned int dim)
      : TestFunction(dim, "ModifiedRidge",
                     std::pair<double, double>{-100.0, 100.0},
                     std::vector<double>(dim, 0.0)) {}

  double value(const std::vector<double> &x) const override {
    if (x.empty())
      return 0.0;

    double term1 = std::fabs(x[0]); // |x1|
    double sum_sq = 0.0;

    for (unsigned int i = 1; i < x.size(); ++i)
      sum_sq += x[i] * x[i];

    return term1 + 2.0 * std::pow(sum_sq, 0.1);
  }
};
/**
 * @class Zakharov
 * @brief Zakharov: f(x) = Σ x_i^2 + (Σ 0.5·i·x_i)^2 + (Σ 0.5·i·x_i)^4.
 * @note Non-separable through the linear coupling terms; minimum at x* = 0.
 */
class Zakharov : public TestFunction {
public:
  Zakharov(unsigned int dim)
      : TestFunction(dim, "Zakharov", std::pair<double, double>{-10.0, 10.0},
                     std::vector<double>(dim, 0.0)) {}

  double value(const std::vector<double> &x) const override {
    if (x.empty())
      return 0.0;

    double sum_sq = 0.0;
    double sum_linear = 0.0;

    for (unsigned int i = 0; i < x.size(); ++i) {
      sum_sq += x[i] * x[i];
      sum_linear += 0.5 * (i + 1) * x[i];
    }

    return sum_sq + std::pow(sum_linear, 2) + std::pow(sum_linear, 4);
  }
};
/**
 * @class ModifiedXinSheYang3
 * @brief Modified Xin-She Yang #3 (scaled): f(x) = 1e4 · (1 + e^{−Σ
 * (x_i/15)^{10}} − 2e^{−Σ x_i^2}) · Π cos^2(x_i).
 * @note Highly multimodal due to product of cos^2 and exponentials.
 */
class ModifiedXinSheYang3 : public TestFunction {
public:
  ModifiedXinSheYang3(unsigned int dim)
      : TestFunction(dim, "ModifiedXin-SheYang",
                     std::pair<double, double>{-20.0, 20.0},
                     std::vector<double>(dim, 0.0)) {}

  double value(const std::vector<double> &x) const override {
    if (x.empty())
      return 0.0;

    double sum_pow10 = 0.0;
    double sum_sq = 0.0;
    double prod_cos2 = 1.0;

    for (double xi : x) {
      sum_pow10 += std::pow(xi / 15.0, 10.0);
      sum_sq += xi * xi;
      prod_cos2 *= std::pow(std::cos(xi), 2.0);
    }

    double inner = std::exp(-sum_pow10) - 2.0 * std::exp(-sum_sq);
    return 1.0e4 * (1.0 + inner) * prod_cos2;
  }
};

/**
 * @class ModifiedXinSheYang5
 * @brief Modified Xin-She Yang #5: f(x) = 1e4·(1 + (Σ sin^2 x_i − e^{Σ
 * x_i^2})·e^{−Σ sin^2(√|x_i|)}).
 * @note Mixes sinusoidal and exponential terms; many local minima.
 */
class ModifiedXinSheYang5 : public TestFunction {
public:
  ModifiedXinSheYang5(unsigned int dim)
      : TestFunction(dim, "ModifiedXin-SheYang5",
                     std::pair<double, double>{-100.0, 100.0},
                     std::vector<double>(dim, 0.0)) {}

  double value(const std::vector<double> &x) const override {
    if (x.empty())
      return 0.0;

    double sum_sin2 = 0.0;         // Σ sin^2(x_i)
    double sum_sq = 0.0;           // Σ x_i^2
    double sum_sin2_sqrtabs = 0.0; // Σ sin^2(√|x_i|)

    for (double xi : x) {
      sum_sin2 += std::pow(std::sin(xi), 2.0);
      sum_sq += xi * xi;
      sum_sin2_sqrtabs += std::pow(std::sin(std::sqrt(std::fabs(xi))), 2.0);
    }

    double inner = sum_sin2 - std::exp(sum_sq);
    return 1.0e4 * (1.0 + inner * std::exp(-sum_sin2_sqrtabs));
  }
};
