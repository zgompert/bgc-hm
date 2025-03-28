// Generated by rstantools.  Do not edit by hand.

/*
    bgchm is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    bgchm is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with bgchm.  If not, see <http://www.gnu.org/licenses/>.
*/
#ifndef MODELS_HPP
#define MODELS_HPP
#define STAN__SERVICES__COMMAND_HPP
#ifndef USE_STANC3
#define USE_STANC3
#endif
#include <rstan/rstaninc.hpp>
// Code generated by stanc v2.32.2
#include <stan/model/model_header.hpp>
namespace model_gencline_one_gl_mix_namespace {
using stan::model::model_base_crtp;
using namespace stan::math;
stan::math::profile_map profiles__;
static constexpr std::array<const char*, 46> locations_array__ =
  {" (found before start of program)",
  " (in 'gencline_one_gl_mix', line 40, column 1 to column 38)",
  " (in 'gencline_one_gl_mix', line 41, column 1 to column 28)",
  " (in 'gencline_one_gl_mix', line 44, column 1 to column 31)",
  " (in 'gencline_one_gl_mix', line 45, column 1 to column 23)",
  " (in 'gencline_one_gl_mix', line 50, column 2 to column 76)",
  " (in 'gencline_one_gl_mix', line 48, column 14 to line 51, column 2)",
  " (in 'gencline_one_gl_mix', line 48, column 1 to line 51, column 2)",
  " (in 'gencline_one_gl_mix', line 53, column 1 to column 46)",
  " (in 'gencline_one_gl_mix', line 54, column 1 to column 41)",
  " (in 'gencline_one_gl_mix', line 27, column 1 to column 7)",
  " (in 'gencline_one_gl_mix', line 28, column 1 to column 7)",
  " (in 'gencline_one_gl_mix', line 29, column 7 to column 8)",
  " (in 'gencline_one_gl_mix', line 29, column 1 to column 37)",
  " (in 'gencline_one_gl_mix', line 30, column 7 to column 8)",
  " (in 'gencline_one_gl_mix', line 30, column 1 to column 37)",
  " (in 'gencline_one_gl_mix', line 31, column 7 to column 8)",
  " (in 'gencline_one_gl_mix', line 31, column 1 to column 37)",
  " (in 'gencline_one_gl_mix', line 32, column 7 to column 8)",
  " (in 'gencline_one_gl_mix', line 32, column 1 to column 40)",
  " (in 'gencline_one_gl_mix', line 33, column 26 to column 27)",
  " (in 'gencline_one_gl_mix', line 33, column 1 to column 31)",
  " (in 'gencline_one_gl_mix', line 34, column 1 to column 27)",
  " (in 'gencline_one_gl_mix', line 35, column 1 to column 27)",
  " (in 'gencline_one_gl_mix', line 36, column 1 to column 18)",
  " (in 'gencline_one_gl_mix', line 37, column 1 to column 18)",
  " (in 'gencline_one_gl_mix', line 3, column 2 to column 11)",
  " (in 'gencline_one_gl_mix', line 4, column 2 to column 43)",
  " (in 'gencline_one_gl_mix', line 5, column 2 to column 13)",
  " (in 'gencline_one_gl_mix', line 2, column 40 to line 6, column 2)",
  " (in 'gencline_one_gl_mix', line 9, column 2 to column 12)",
  " (in 'gencline_one_gl_mix', line 10, column 2 to column 11)",
  " (in 'gencline_one_gl_mix', line 11, column 2 to column 28)",
  " (in 'gencline_one_gl_mix', line 21, column 3 to column 12)",
  " (in 'gencline_one_gl_mix', line 20, column 8 to line 22, column 3)",
  " (in 'gencline_one_gl_mix', line 17, column 3 to column 58)",
  " (in 'gencline_one_gl_mix', line 19, column 3 to column 69)",
  " (in 'gencline_one_gl_mix', line 16, column 19 to line 20, column 3)",
  " (in 'gencline_one_gl_mix', line 16, column 9 to line 22, column 3)",
  " (in 'gencline_one_gl_mix', line 13, column 14 to column 108)",
  " (in 'gencline_one_gl_mix', line 14, column 14 to column 128)",
  " (in 'gencline_one_gl_mix', line 15, column 14 to column 111)",
  " (in 'gencline_one_gl_mix', line 12, column 11 to line 16, column 3)",
  " (in 'gencline_one_gl_mix', line 12, column 2 to line 22, column 3)",
  " (in 'gencline_one_gl_mix', line 23, column 2 to column 14)",
  " (in 'gencline_one_gl_mix', line 8, column 97 to line 24, column 2)"};
template <typename T0__, typename T1__, typename T2__,
          stan::require_all_t<stan::is_stan_scalar<T0__>,
                              stan::is_stan_scalar<T1__>,
                              stan::is_stan_scalar<T2__>>* = nullptr>
stan::promote_args_t<T0__, T1__, T2__>
calc_phi(const T0__& h, const T1__& vv, const T2__& uu, std::ostream*
         pstream__);
template <typename T0__, typename T1__, typename T2__, typename T3__,
          typename T4__, typename T5__, typename T6__, typename T7__,
          typename T8__,
          stan::require_all_t<stan::is_stan_scalar<T0__>,
                              stan::is_stan_scalar<T1__>,
                              stan::is_stan_scalar<T2__>,
                              stan::is_stan_scalar<T3__>,
                              stan::is_stan_scalar<T4__>,
                              stan::is_stan_scalar<T5__>,
                              stan::is_stan_scalar<T6__>,
                              stan::is_stan_scalar<T7__>,
                              stan::is_stan_scalar<T8__>>* = nullptr>
stan::promote_args_t<T0__, T1__, T2__, T3__, T4__,
  stan::promote_args_t<T5__, T6__, T7__, T8__>>
calc_lik(const T0__& gl0, const T1__& gl1, const T2__& gl2, const T3__& p0,
         const T4__& p1, const T5__& h, const T6__& vv, const T7__& uu,
         const T8__& pl, std::ostream* pstream__);
template <typename T0__, typename T1__, typename T2__,
          stan::require_all_t<stan::is_stan_scalar<T0__>,
                              stan::is_stan_scalar<T1__>,
                              stan::is_stan_scalar<T2__>>*>
stan::promote_args_t<T0__, T1__, T2__>
calc_phi(const T0__& h, const T1__& vv, const T2__& uu, std::ostream*
         pstream__) {
  using local_scalar_t__ = stan::promote_args_t<T0__, T1__, T2__>;
  int current_statement__ = 0;
  static constexpr bool propto__ = true;
  // suppress unused var warning
  (void) propto__;
  local_scalar_t__ DUMMY_VAR__(std::numeric_limits<double>::quiet_NaN());
  // suppress unused var warning
  (void) DUMMY_VAR__;
  try {
    local_scalar_t__ phi = DUMMY_VAR__;
    current_statement__ = 27;
    phi = (stan::math::pow(h, vv) / (stan::math::pow(h, vv) +
      (stan::math::pow((1 - h), vv) * stan::math::exp(uu))));
    current_statement__ = 28;
    return phi;
  } catch (const std::exception& e) {
    stan::lang::rethrow_located(e, locations_array__[current_statement__]);
  }
}
template <typename T0__, typename T1__, typename T2__, typename T3__,
          typename T4__, typename T5__, typename T6__, typename T7__,
          typename T8__,
          stan::require_all_t<stan::is_stan_scalar<T0__>,
                              stan::is_stan_scalar<T1__>,
                              stan::is_stan_scalar<T2__>,
                              stan::is_stan_scalar<T3__>,
                              stan::is_stan_scalar<T4__>,
                              stan::is_stan_scalar<T5__>,
                              stan::is_stan_scalar<T6__>,
                              stan::is_stan_scalar<T7__>,
                              stan::is_stan_scalar<T8__>>*>
stan::promote_args_t<T0__, T1__, T2__, T3__, T4__,
  stan::promote_args_t<T5__, T6__, T7__, T8__>>
calc_lik(const T0__& gl0, const T1__& gl1, const T2__& gl2, const T3__& p0,
         const T4__& p1, const T5__& h, const T6__& vv, const T7__& uu,
         const T8__& pl, std::ostream* pstream__) {
  using local_scalar_t__ = stan::promote_args_t<T0__, T1__, T2__, T3__, T4__,
                             stan::promote_args_t<T5__, T6__, T7__, T8__>>;
  int current_statement__ = 0;
  static constexpr bool propto__ = true;
  // suppress unused var warning
  (void) propto__;
  local_scalar_t__ DUMMY_VAR__(std::numeric_limits<double>::quiet_NaN());
  // suppress unused var warning
  (void) DUMMY_VAR__;
  try {
    local_scalar_t__ prob = DUMMY_VAR__;
    local_scalar_t__ phi = DUMMY_VAR__;
    current_statement__ = 32;
    phi = calc_phi(h, vv, uu, pstream__);
    current_statement__ = 43;
    if (stan::math::logical_eq(pl, 2)) {
      current_statement__ = 39;
      prob = ((stan::math::log(gl0) +
        stan::math::log(((phi * (1 - p1)) + ((1 - phi) * (1 - p0))))) +
        stan::math::log(((phi * (1 - p1)) + ((1 - phi) * (1 - p0)))));
      current_statement__ = 40;
      prob = stan::math::log_sum_exp(prob, (((stan::math::log(gl1) +
               stan::math::log(2)) +
               stan::math::log(((phi * (1 - p1)) + ((1 - phi) * (1 - p0)))))
               + stan::math::log(((phi * p1) + ((1 - phi) * p0)))));
      current_statement__ = 41;
      prob = stan::math::log_sum_exp(prob, ((stan::math::log(gl2) +
               stan::math::log(((phi * p1) + ((1 - phi) * p0)))) +
               stan::math::log(((phi * p1) + ((1 - phi) * p0)))));
    } else {
      current_statement__ = 38;
      if (stan::math::logical_eq(pl, 1)) {
        current_statement__ = 35;
        prob = (stan::math::log(gl0) +
          stan::math::log(((phi * (1 - p1)) + ((1 - phi) * (1 - p0)))));
        current_statement__ = 36;
        prob = stan::math::log_sum_exp(prob, (stan::math::log(gl1) +
                 stan::math::log(((phi * p1) + ((1 - phi) * p0)))));
      } else {
        current_statement__ = 33;
        prob = 0;
      }
    }
    current_statement__ = 44;
    return prob;
  } catch (const std::exception& e) {
    stan::lang::rethrow_located(e, locations_array__[current_statement__]);
  }
}
#include <stan_meta_header.hpp>
class model_gencline_one_gl_mix final : public model_base_crtp<model_gencline_one_gl_mix> {
private:
  int L;
  int N;
  std::vector<double> GL0;
  std::vector<double> GL1;
  std::vector<double> GL2;
  std::vector<double> ploidy;
  Eigen::Matrix<double,-1,1> H_data__;
  double P0;
  double P1;
  double sc;
  double sv;
  Eigen::Map<Eigen::Matrix<double,-1,1>> H{nullptr, 0};
public:
  ~model_gencline_one_gl_mix() {}
  model_gencline_one_gl_mix(stan::io::var_context& context__, unsigned int
                            random_seed__ = 0, std::ostream*
                            pstream__ = nullptr) : model_base_crtp(0) {
    int current_statement__ = 0;
    using local_scalar_t__ = double;
    boost::ecuyer1988 base_rng__ =
      stan::services::util::create_rng(random_seed__, 0);
    // suppress unused var warning
    (void) base_rng__;
    static constexpr const char* function__ =
      "model_gencline_one_gl_mix_namespace::model_gencline_one_gl_mix";
    // suppress unused var warning
    (void) function__;
    local_scalar_t__ DUMMY_VAR__(std::numeric_limits<double>::quiet_NaN());
    // suppress unused var warning
    (void) DUMMY_VAR__;
    try {
      int pos__ = std::numeric_limits<int>::min();
      pos__ = 1;
      current_statement__ = 10;
      context__.validate_dims("data initialization", "L", "int",
        std::vector<size_t>{});
      L = std::numeric_limits<int>::min();
      current_statement__ = 10;
      L = context__.vals_i("L")[(1 - 1)];
      current_statement__ = 11;
      context__.validate_dims("data initialization", "N", "int",
        std::vector<size_t>{});
      N = std::numeric_limits<int>::min();
      current_statement__ = 11;
      N = context__.vals_i("N")[(1 - 1)];
      current_statement__ = 12;
      stan::math::validate_non_negative_index("GL0", "N", N);
      current_statement__ = 13;
      context__.validate_dims("data initialization", "GL0", "double",
        std::vector<size_t>{static_cast<size_t>(N)});
      GL0 = std::vector<double>(N, std::numeric_limits<double>::quiet_NaN());
      current_statement__ = 13;
      GL0 = context__.vals_r("GL0");
      current_statement__ = 13;
      stan::math::check_greater_or_equal(function__, "GL0", GL0, 0);
      current_statement__ = 13;
      stan::math::check_less_or_equal(function__, "GL0", GL0, 2);
      current_statement__ = 14;
      stan::math::validate_non_negative_index("GL1", "N", N);
      current_statement__ = 15;
      context__.validate_dims("data initialization", "GL1", "double",
        std::vector<size_t>{static_cast<size_t>(N)});
      GL1 = std::vector<double>(N, std::numeric_limits<double>::quiet_NaN());
      current_statement__ = 15;
      GL1 = context__.vals_r("GL1");
      current_statement__ = 15;
      stan::math::check_greater_or_equal(function__, "GL1", GL1, 0);
      current_statement__ = 15;
      stan::math::check_less_or_equal(function__, "GL1", GL1, 2);
      current_statement__ = 16;
      stan::math::validate_non_negative_index("GL2", "N", N);
      current_statement__ = 17;
      context__.validate_dims("data initialization", "GL2", "double",
        std::vector<size_t>{static_cast<size_t>(N)});
      GL2 = std::vector<double>(N, std::numeric_limits<double>::quiet_NaN());
      current_statement__ = 17;
      GL2 = context__.vals_r("GL2");
      current_statement__ = 17;
      stan::math::check_greater_or_equal(function__, "GL2", GL2, 0);
      current_statement__ = 17;
      stan::math::check_less_or_equal(function__, "GL2", GL2, 2);
      current_statement__ = 18;
      stan::math::validate_non_negative_index("ploidy", "N", N);
      current_statement__ = 19;
      context__.validate_dims("data initialization", "ploidy", "double",
        std::vector<size_t>{static_cast<size_t>(N)});
      ploidy = std::vector<double>(N,
                 std::numeric_limits<double>::quiet_NaN());
      current_statement__ = 19;
      ploidy = context__.vals_r("ploidy");
      current_statement__ = 19;
      stan::math::check_greater_or_equal(function__, "ploidy", ploidy, 0);
      current_statement__ = 19;
      stan::math::check_less_or_equal(function__, "ploidy", ploidy, 2);
      current_statement__ = 20;
      stan::math::validate_non_negative_index("H", "N", N);
      current_statement__ = 21;
      context__.validate_dims("data initialization", "H", "double",
        std::vector<size_t>{static_cast<size_t>(N)});
      H_data__ = Eigen::Matrix<double,-1,1>::Constant(N,
                   std::numeric_limits<double>::quiet_NaN());
      new (&H) Eigen::Map<Eigen::Matrix<double,-1,1>>(H_data__.data(), N);
      {
        std::vector<local_scalar_t__> H_flat__;
        current_statement__ = 21;
        H_flat__ = context__.vals_r("H");
        current_statement__ = 21;
        pos__ = 1;
        current_statement__ = 21;
        for (int sym1__ = 1; sym1__ <= N; ++sym1__) {
          current_statement__ = 21;
          stan::model::assign(H, H_flat__[(pos__ - 1)],
            "assigning variable H", stan::model::index_uni(sym1__));
          current_statement__ = 21;
          pos__ = (pos__ + 1);
        }
      }
      current_statement__ = 21;
      stan::math::check_greater_or_equal(function__, "H", H, 0);
      current_statement__ = 21;
      stan::math::check_less_or_equal(function__, "H", H, 1);
      current_statement__ = 22;
      context__.validate_dims("data initialization", "P0", "double",
        std::vector<size_t>{});
      P0 = std::numeric_limits<double>::quiet_NaN();
      current_statement__ = 22;
      P0 = context__.vals_r("P0")[(1 - 1)];
      current_statement__ = 22;
      stan::math::check_greater_or_equal(function__, "P0", P0, 0);
      current_statement__ = 22;
      stan::math::check_less_or_equal(function__, "P0", P0, 1);
      current_statement__ = 23;
      context__.validate_dims("data initialization", "P1", "double",
        std::vector<size_t>{});
      P1 = std::numeric_limits<double>::quiet_NaN();
      current_statement__ = 23;
      P1 = context__.vals_r("P1")[(1 - 1)];
      current_statement__ = 23;
      stan::math::check_greater_or_equal(function__, "P1", P1, 0);
      current_statement__ = 23;
      stan::math::check_less_or_equal(function__, "P1", P1, 1);
      current_statement__ = 24;
      context__.validate_dims("data initialization", "sc", "double",
        std::vector<size_t>{});
      sc = std::numeric_limits<double>::quiet_NaN();
      current_statement__ = 24;
      sc = context__.vals_r("sc")[(1 - 1)];
      current_statement__ = 24;
      stan::math::check_greater_or_equal(function__, "sc", sc, 0);
      current_statement__ = 25;
      context__.validate_dims("data initialization", "sv", "double",
        std::vector<size_t>{});
      sv = std::numeric_limits<double>::quiet_NaN();
      current_statement__ = 25;
      sv = context__.vals_r("sv")[(1 - 1)];
      current_statement__ = 25;
      stan::math::check_greater_or_equal(function__, "sv", sv, 0);
    } catch (const std::exception& e) {
      stan::lang::rethrow_located(e, locations_array__[current_statement__]);
    }
    num_params_r__ = 1 + 1;
  }
  inline std::string model_name() const final {
    return "model_gencline_one_gl_mix";
  }
  inline std::vector<std::string> model_compile_info() const noexcept {
    return std::vector<std::string>{"stanc_version = stanc3 v2.32.2",
             "stancflags = --allow-undefined"};
  }
  template <bool propto__, bool jacobian__, typename VecR, typename VecI,
            stan::require_vector_like_t<VecR>* = nullptr,
            stan::require_vector_like_vt<std::is_integral, VecI>* = nullptr>
  inline stan::scalar_type_t<VecR>
  log_prob_impl(VecR& params_r__, VecI& params_i__, std::ostream*
                pstream__ = nullptr) const {
    using T__ = stan::scalar_type_t<VecR>;
    using local_scalar_t__ = T__;
    T__ lp__(0.0);
    stan::math::accumulator<T__> lp_accum__;
    stan::io::deserializer<local_scalar_t__> in__(params_r__, params_i__);
    int current_statement__ = 0;
    local_scalar_t__ DUMMY_VAR__(std::numeric_limits<double>::quiet_NaN());
    // suppress unused var warning
    (void) DUMMY_VAR__;
    static constexpr const char* function__ =
      "model_gencline_one_gl_mix_namespace::log_prob";
    // suppress unused var warning
    (void) function__;
    try {
      local_scalar_t__ center = DUMMY_VAR__;
      current_statement__ = 1;
      center = in__.template read_constrain_lub<local_scalar_t__,
                 jacobian__>(0.001, 0.999, lp__);
      local_scalar_t__ v = DUMMY_VAR__;
      current_statement__ = 2;
      v = in__.template read_constrain_lub<local_scalar_t__, jacobian__>(0.1,
            10, lp__);
      local_scalar_t__ u = DUMMY_VAR__;
      current_statement__ = 4;
      u = (stan::math::logit(center) * v);
      current_statement__ = 3;
      stan::math::check_greater_or_equal(function__, "u", u, -100);
      current_statement__ = 3;
      stan::math::check_less_or_equal(function__, "u", u, 100);
      {
        current_statement__ = 7;
        for (int j = 1; j <= N; ++j) {
          current_statement__ = 5;
          lp_accum__.add(calc_lik(
                           stan::model::rvalue(GL0, "GL0",
                             stan::model::index_uni(j)),
                           stan::model::rvalue(GL1, "GL1",
                             stan::model::index_uni(j)),
                           stan::model::rvalue(GL2, "GL2",
                             stan::model::index_uni(j)), P0, P1,
                           stan::model::rvalue(H, "H",
                             stan::model::index_uni(j)), v, u,
                           stan::model::rvalue(ploidy, "ploidy",
                             stan::model::index_uni(j)), pstream__));
        }
        current_statement__ = 8;
        lp_accum__.add(stan::math::normal_lpdf<false>(
                         stan::math::logit(center), 0, sc));
        current_statement__ = 9;
        lp_accum__.add(stan::math::normal_lpdf<false>(stan::math::log10(v),
                         0, sv));
      }
    } catch (const std::exception& e) {
      stan::lang::rethrow_located(e, locations_array__[current_statement__]);
    }
    lp_accum__.add(lp__);
    return lp_accum__.sum();
  }
  template <typename RNG, typename VecR, typename VecI, typename VecVar,
            stan::require_vector_like_vt<std::is_floating_point,
            VecR>* = nullptr, stan::require_vector_like_vt<std::is_integral,
            VecI>* = nullptr, stan::require_vector_vt<std::is_floating_point,
            VecVar>* = nullptr>
  inline void
  write_array_impl(RNG& base_rng__, VecR& params_r__, VecI& params_i__,
                   VecVar& vars__, const bool
                   emit_transformed_parameters__ = true, const bool
                   emit_generated_quantities__ = true, std::ostream*
                   pstream__ = nullptr) const {
    using local_scalar_t__ = double;
    stan::io::deserializer<local_scalar_t__> in__(params_r__, params_i__);
    stan::io::serializer<local_scalar_t__> out__(vars__);
    static constexpr bool propto__ = true;
    // suppress unused var warning
    (void) propto__;
    double lp__ = 0.0;
    // suppress unused var warning
    (void) lp__;
    int current_statement__ = 0;
    stan::math::accumulator<double> lp_accum__;
    local_scalar_t__ DUMMY_VAR__(std::numeric_limits<double>::quiet_NaN());
    // suppress unused var warning
    (void) DUMMY_VAR__;
    constexpr bool jacobian__ = false;
    static constexpr const char* function__ =
      "model_gencline_one_gl_mix_namespace::write_array";
    // suppress unused var warning
    (void) function__;
    try {
      double center = std::numeric_limits<double>::quiet_NaN();
      current_statement__ = 1;
      center = in__.template read_constrain_lub<local_scalar_t__,
                 jacobian__>(0.001, 0.999, lp__);
      double v = std::numeric_limits<double>::quiet_NaN();
      current_statement__ = 2;
      v = in__.template read_constrain_lub<local_scalar_t__, jacobian__>(0.1,
            10, lp__);
      double u = std::numeric_limits<double>::quiet_NaN();
      out__.write(center);
      out__.write(v);
      if (stan::math::logical_negation(
            (stan::math::primitive_value(emit_transformed_parameters__) ||
            stan::math::primitive_value(emit_generated_quantities__)))) {
        return ;
      }
      current_statement__ = 4;
      u = (stan::math::logit(center) * v);
      current_statement__ = 3;
      stan::math::check_greater_or_equal(function__, "u", u, -100);
      current_statement__ = 3;
      stan::math::check_less_or_equal(function__, "u", u, 100);
      if (emit_transformed_parameters__) {
        out__.write(u);
      }
      if (stan::math::logical_negation(emit_generated_quantities__)) {
        return ;
      }
    } catch (const std::exception& e) {
      stan::lang::rethrow_located(e, locations_array__[current_statement__]);
    }
  }
  template <typename VecVar, typename VecI,
            stan::require_vector_t<VecVar>* = nullptr,
            stan::require_vector_like_vt<std::is_integral, VecI>* = nullptr>
  inline void
  unconstrain_array_impl(const VecVar& params_r__, const VecI& params_i__,
                         VecVar& vars__, std::ostream* pstream__ = nullptr) const {
    using local_scalar_t__ = double;
    stan::io::deserializer<local_scalar_t__> in__(params_r__, params_i__);
    stan::io::serializer<local_scalar_t__> out__(vars__);
    int current_statement__ = 0;
    local_scalar_t__ DUMMY_VAR__(std::numeric_limits<double>::quiet_NaN());
    // suppress unused var warning
    (void) DUMMY_VAR__;
    try {
      int pos__ = std::numeric_limits<int>::min();
      pos__ = 1;
      local_scalar_t__ center = DUMMY_VAR__;
      current_statement__ = 1;
      center = in__.read<local_scalar_t__>();
      out__.write_free_lub(0.001, 0.999, center);
      local_scalar_t__ v = DUMMY_VAR__;
      current_statement__ = 2;
      v = in__.read<local_scalar_t__>();
      out__.write_free_lub(0.1, 10, v);
    } catch (const std::exception& e) {
      stan::lang::rethrow_located(e, locations_array__[current_statement__]);
    }
  }
  template <typename VecVar, stan::require_vector_t<VecVar>* = nullptr>
  inline void
  transform_inits_impl(const stan::io::var_context& context__, VecVar&
                       vars__, std::ostream* pstream__ = nullptr) const {
    using local_scalar_t__ = double;
    stan::io::serializer<local_scalar_t__> out__(vars__);
    int current_statement__ = 0;
    local_scalar_t__ DUMMY_VAR__(std::numeric_limits<double>::quiet_NaN());
    // suppress unused var warning
    (void) DUMMY_VAR__;
    try {
      current_statement__ = 1;
      context__.validate_dims("parameter initialization", "center", "double",
        std::vector<size_t>{});
      current_statement__ = 2;
      context__.validate_dims("parameter initialization", "v", "double",
        std::vector<size_t>{});
      int pos__ = std::numeric_limits<int>::min();
      pos__ = 1;
      local_scalar_t__ center = DUMMY_VAR__;
      current_statement__ = 1;
      center = context__.vals_r("center")[(1 - 1)];
      out__.write_free_lub(0.001, 0.999, center);
      local_scalar_t__ v = DUMMY_VAR__;
      current_statement__ = 2;
      v = context__.vals_r("v")[(1 - 1)];
      out__.write_free_lub(0.1, 10, v);
    } catch (const std::exception& e) {
      stan::lang::rethrow_located(e, locations_array__[current_statement__]);
    }
  }
  inline void
  get_param_names(std::vector<std::string>& names__, const bool
                  emit_transformed_parameters__ = true, const bool
                  emit_generated_quantities__ = true) const {
    names__ = std::vector<std::string>{"center", "v"};
    if (emit_transformed_parameters__) {
      std::vector<std::string> temp{"u"};
      names__.reserve(names__.size() + temp.size());
      names__.insert(names__.end(), temp.begin(), temp.end());
    }
    if (emit_generated_quantities__) {}
  }
  inline void
  get_dims(std::vector<std::vector<size_t>>& dimss__, const bool
           emit_transformed_parameters__ = true, const bool
           emit_generated_quantities__ = true) const {
    dimss__ = std::vector<std::vector<size_t>>{std::vector<size_t>{},
                std::vector<size_t>{}};
    if (emit_transformed_parameters__) {
      std::vector<std::vector<size_t>> temp{std::vector<size_t>{}};
      dimss__.reserve(dimss__.size() + temp.size());
      dimss__.insert(dimss__.end(), temp.begin(), temp.end());
    }
    if (emit_generated_quantities__) {}
  }
  inline void
  constrained_param_names(std::vector<std::string>& param_names__, bool
                          emit_transformed_parameters__ = true, bool
                          emit_generated_quantities__ = true) const final {
    param_names__.emplace_back(std::string() + "center");
    param_names__.emplace_back(std::string() + "v");
    if (emit_transformed_parameters__) {
      param_names__.emplace_back(std::string() + "u");
    }
    if (emit_generated_quantities__) {}
  }
  inline void
  unconstrained_param_names(std::vector<std::string>& param_names__, bool
                            emit_transformed_parameters__ = true, bool
                            emit_generated_quantities__ = true) const final {
    param_names__.emplace_back(std::string() + "center");
    param_names__.emplace_back(std::string() + "v");
    if (emit_transformed_parameters__) {
      param_names__.emplace_back(std::string() + "u");
    }
    if (emit_generated_quantities__) {}
  }
  inline std::string get_constrained_sizedtypes() const {
    return std::string("[{\"name\":\"center\",\"type\":{\"name\":\"real\"},\"block\":\"parameters\"},{\"name\":\"v\",\"type\":{\"name\":\"real\"},\"block\":\"parameters\"},{\"name\":\"u\",\"type\":{\"name\":\"real\"},\"block\":\"transformed_parameters\"}]");
  }
  inline std::string get_unconstrained_sizedtypes() const {
    return std::string("[{\"name\":\"center\",\"type\":{\"name\":\"real\"},\"block\":\"parameters\"},{\"name\":\"v\",\"type\":{\"name\":\"real\"},\"block\":\"parameters\"},{\"name\":\"u\",\"type\":{\"name\":\"real\"},\"block\":\"transformed_parameters\"}]");
  }
  // Begin method overload boilerplate
  template <typename RNG> inline void
  write_array(RNG& base_rng, Eigen::Matrix<double,-1,1>& params_r,
              Eigen::Matrix<double,-1,1>& vars, const bool
              emit_transformed_parameters = true, const bool
              emit_generated_quantities = true, std::ostream*
              pstream = nullptr) const {
    const size_t num_params__ = (1 + 1);
    const size_t num_transformed = emit_transformed_parameters * (1);
    const size_t num_gen_quantities = emit_generated_quantities * (0);
    const size_t num_to_write = num_params__ + num_transformed +
      num_gen_quantities;
    std::vector<int> params_i;
    vars = Eigen::Matrix<double,-1,1>::Constant(num_to_write,
             std::numeric_limits<double>::quiet_NaN());
    write_array_impl(base_rng, params_r, params_i, vars,
      emit_transformed_parameters, emit_generated_quantities, pstream);
  }
  template <typename RNG> inline void
  write_array(RNG& base_rng, std::vector<double>& params_r, std::vector<int>&
              params_i, std::vector<double>& vars, bool
              emit_transformed_parameters = true, bool
              emit_generated_quantities = true, std::ostream*
              pstream = nullptr) const {
    const size_t num_params__ = (1 + 1);
    const size_t num_transformed = emit_transformed_parameters * (1);
    const size_t num_gen_quantities = emit_generated_quantities * (0);
    const size_t num_to_write = num_params__ + num_transformed +
      num_gen_quantities;
    vars = std::vector<double>(num_to_write,
             std::numeric_limits<double>::quiet_NaN());
    write_array_impl(base_rng, params_r, params_i, vars,
      emit_transformed_parameters, emit_generated_quantities, pstream);
  }
  template <bool propto__, bool jacobian__, typename T_> inline T_
  log_prob(Eigen::Matrix<T_,-1,1>& params_r, std::ostream* pstream = nullptr) const {
    Eigen::Matrix<int,-1,1> params_i;
    return log_prob_impl<propto__, jacobian__>(params_r, params_i, pstream);
  }
  template <bool propto__, bool jacobian__, typename T_> inline T_
  log_prob(std::vector<T_>& params_r, std::vector<int>& params_i,
           std::ostream* pstream = nullptr) const {
    return log_prob_impl<propto__, jacobian__>(params_r, params_i, pstream);
  }
  inline void
  transform_inits(const stan::io::var_context& context,
                  Eigen::Matrix<double,-1,1>& params_r, std::ostream*
                  pstream = nullptr) const final {
    std::vector<double> params_r_vec(params_r.size());
    std::vector<int> params_i;
    transform_inits(context, params_i, params_r_vec, pstream);
    params_r = Eigen::Map<Eigen::Matrix<double,-1,1>>(params_r_vec.data(),
                 params_r_vec.size());
  }
  inline void
  transform_inits(const stan::io::var_context& context, std::vector<int>&
                  params_i, std::vector<double>& vars, std::ostream*
                  pstream__ = nullptr) const {
    vars.resize(num_params_r__);
    transform_inits_impl(context, vars, pstream__);
  }
  inline void
  unconstrain_array(const std::vector<double>& params_constrained,
                    std::vector<double>& params_unconstrained, std::ostream*
                    pstream = nullptr) const {
    const std::vector<int> params_i;
    params_unconstrained = std::vector<double>(num_params_r__,
                             std::numeric_limits<double>::quiet_NaN());
    unconstrain_array_impl(params_constrained, params_i,
      params_unconstrained, pstream);
  }
  inline void
  unconstrain_array(const Eigen::Matrix<double,-1,1>& params_constrained,
                    Eigen::Matrix<double,-1,1>& params_unconstrained,
                    std::ostream* pstream = nullptr) const {
    const std::vector<int> params_i;
    params_unconstrained = Eigen::Matrix<double,-1,1>::Constant(num_params_r__,
                             std::numeric_limits<double>::quiet_NaN());
    unconstrain_array_impl(params_constrained, params_i,
      params_unconstrained, pstream);
  }
};
}
using stan_model = model_gencline_one_gl_mix_namespace::model_gencline_one_gl_mix;
#ifndef USING_R
// Boilerplate
stan::model::model_base&
new_model(stan::io::var_context& data_context, unsigned int seed,
          std::ostream* msg_stream) {
  stan_model* m = new stan_model(data_context, seed, msg_stream);
  return *m;
}
stan::math::profile_map& get_stan_profile_data() {
  return model_gencline_one_gl_mix_namespace::profiles__;
}
#endif
#endif
