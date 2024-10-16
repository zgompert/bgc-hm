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
#include <rstan/rstaninc.hpp>
// Code generated by Stan version 2.21.0
#include <stan/model/model_header.hpp>
namespace model_popp_gl_mix_namespace {
using std::istream;
using std::string;
using std::stringstream;
using std::vector;
using stan::io::dump;
using stan::math::lgamma;
using stan::model::prob_grad;
using namespace stan::math;
static int current_statement_begin__;
stan::io::program_reader prog_reader__() {
    stan::io::program_reader reader;
    reader.add_event(0, 0, "start", "model_popp_gl_mix");
    reader.add_event(49, 47, "end", "model_popp_gl_mix");
    return reader;
}
template <typename T0__, typename T1__, typename T2__, typename T3__, typename T4__>
typename boost::math::tools::promote_args<T0__, T1__, T2__, T3__, typename boost::math::tools::promote_args<T4__>::type>::type
calc_lik(const T0__& gl0,
             const T1__& gl1,
             const T2__& gl2,
             const T3__& p,
             const T4__& pl, std::ostream* pstream__) {
    typedef typename boost::math::tools::promote_args<T0__, T1__, T2__, T3__, typename boost::math::tools::promote_args<T4__>::type>::type local_scalar_t__;
    typedef local_scalar_t__ fun_return_scalar_t__;
    const static bool propto__ = true;
    (void) propto__;
        local_scalar_t__ DUMMY_VAR__(std::numeric_limits<double>::quiet_NaN());
        (void) DUMMY_VAR__;  // suppress unused var warning
    int current_statement_begin__ = -1;
    try {
        {
        current_statement_begin__ = 4;
        local_scalar_t__ prob(DUMMY_VAR__);
        (void) prob;  // dummy to suppress unused var warning
        stan::math::initialize(prob, DUMMY_VAR__);
        stan::math::fill(prob, DUMMY_VAR__);
        current_statement_begin__ = 5;
        if (as_bool(logical_eq(pl, 2))) {
            current_statement_begin__ = 6;
            stan::math::assign(prob, ((stan::math::log(gl0) + stan::math::log((1 - p))) + stan::math::log((1 - p))));
            current_statement_begin__ = 7;
            stan::math::assign(prob, log_sum_exp(prob, ((((stan::math::log(gl1) + stan::math::log(2)) + stan::math::log(p)) + stan::math::log((1 - p))) + stan::math::log(2))));
            current_statement_begin__ = 8;
            stan::math::assign(prob, log_sum_exp(prob, ((stan::math::log(gl2) + stan::math::log(p)) + stan::math::log(p))));
        } else if (as_bool(logical_eq(pl, 1))) {
            current_statement_begin__ = 10;
            stan::math::assign(prob, (stan::math::log(gl0) + stan::math::log((1 - p))));
            current_statement_begin__ = 11;
            stan::math::assign(prob, log_sum_exp(prob, (stan::math::log(gl1) + stan::math::log(p))));
        } else {
            current_statement_begin__ = 13;
            stan::math::assign(prob, 0);
        }
        current_statement_begin__ = 15;
        return stan::math::promote_scalar<fun_return_scalar_t__>(prob);
        }
    } catch (const std::exception& e) {
        stan::lang::rethrow_located(e, current_statement_begin__, prog_reader__());
        // Next line prevents compiler griping about no return
        throw std::runtime_error("*** IF YOU SEE THIS, PLEASE REPORT A BUG ***");
    }
}
struct calc_lik_functor__ {
    template <typename T0__, typename T1__, typename T2__, typename T3__, typename T4__>
        typename boost::math::tools::promote_args<T0__, T1__, T2__, T3__, typename boost::math::tools::promote_args<T4__>::type>::type
    operator()(const T0__& gl0,
             const T1__& gl1,
             const T2__& gl2,
             const T3__& p,
             const T4__& pl, std::ostream* pstream__) const {
        return calc_lik(gl0, gl1, gl2, p, pl, pstream__);
    }
};
#include <stan_meta_header.hpp>
class model_popp_gl_mix
  : public stan::model::model_base_crtp<model_popp_gl_mix> {
private:
        int L;
        int N;
        int J;
        std::vector<std::vector<double> > GL0;
        std::vector<std::vector<double> > GL1;
        std::vector<std::vector<double> > GL2;
        std::vector<double> ploidy;
        std::vector<int> pids;
public:
    model_popp_gl_mix(stan::io::var_context& context__,
        std::ostream* pstream__ = 0)
        : model_base_crtp(0) {
        ctor_body(context__, 0, pstream__);
    }
    model_popp_gl_mix(stan::io::var_context& context__,
        unsigned int random_seed__,
        std::ostream* pstream__ = 0)
        : model_base_crtp(0) {
        ctor_body(context__, random_seed__, pstream__);
    }
    void ctor_body(stan::io::var_context& context__,
                   unsigned int random_seed__,
                   std::ostream* pstream__) {
        typedef double local_scalar_t__;
        boost::ecuyer1988 base_rng__ =
          stan::services::util::create_rng(random_seed__, 0);
        (void) base_rng__;  // suppress unused var warning
        current_statement_begin__ = -1;
        static const char* function__ = "model_popp_gl_mix_namespace::model_popp_gl_mix";
        (void) function__;  // dummy to suppress unused var warning
        size_t pos__;
        (void) pos__;  // dummy to suppress unused var warning
        std::vector<int> vals_i__;
        std::vector<double> vals_r__;
        local_scalar_t__ DUMMY_VAR__(std::numeric_limits<double>::quiet_NaN());
        (void) DUMMY_VAR__;  // suppress unused var warning
        try {
            // initialize data block variables from context__
            current_statement_begin__ = 19;
            context__.validate_dims("data initialization", "L", "int", context__.to_vec());
            L = int(0);
            vals_i__ = context__.vals_i("L");
            pos__ = 0;
            L = vals_i__[pos__++];
            current_statement_begin__ = 20;
            context__.validate_dims("data initialization", "N", "int", context__.to_vec());
            N = int(0);
            vals_i__ = context__.vals_i("N");
            pos__ = 0;
            N = vals_i__[pos__++];
            current_statement_begin__ = 21;
            context__.validate_dims("data initialization", "J", "int", context__.to_vec());
            J = int(0);
            vals_i__ = context__.vals_i("J");
            pos__ = 0;
            J = vals_i__[pos__++];
            current_statement_begin__ = 22;
            validate_non_negative_index("GL0", "N", N);
            validate_non_negative_index("GL0", "L", L);
            context__.validate_dims("data initialization", "GL0", "double", context__.to_vec(N,L));
            GL0 = std::vector<std::vector<double> >(N, std::vector<double>(L, double(0)));
            vals_r__ = context__.vals_r("GL0");
            pos__ = 0;
            size_t GL0_k_0_max__ = N;
            size_t GL0_k_1_max__ = L;
            for (size_t k_1__ = 0; k_1__ < GL0_k_1_max__; ++k_1__) {
                for (size_t k_0__ = 0; k_0__ < GL0_k_0_max__; ++k_0__) {
                    GL0[k_0__][k_1__] = vals_r__[pos__++];
                }
            }
            size_t GL0_i_0_max__ = N;
            size_t GL0_i_1_max__ = L;
            for (size_t i_0__ = 0; i_0__ < GL0_i_0_max__; ++i_0__) {
                for (size_t i_1__ = 0; i_1__ < GL0_i_1_max__; ++i_1__) {
                    check_greater_or_equal(function__, "GL0[i_0__][i_1__]", GL0[i_0__][i_1__], 0);
                    check_less_or_equal(function__, "GL0[i_0__][i_1__]", GL0[i_0__][i_1__], 1);
                }
            }
            current_statement_begin__ = 23;
            validate_non_negative_index("GL1", "N", N);
            validate_non_negative_index("GL1", "L", L);
            context__.validate_dims("data initialization", "GL1", "double", context__.to_vec(N,L));
            GL1 = std::vector<std::vector<double> >(N, std::vector<double>(L, double(0)));
            vals_r__ = context__.vals_r("GL1");
            pos__ = 0;
            size_t GL1_k_0_max__ = N;
            size_t GL1_k_1_max__ = L;
            for (size_t k_1__ = 0; k_1__ < GL1_k_1_max__; ++k_1__) {
                for (size_t k_0__ = 0; k_0__ < GL1_k_0_max__; ++k_0__) {
                    GL1[k_0__][k_1__] = vals_r__[pos__++];
                }
            }
            size_t GL1_i_0_max__ = N;
            size_t GL1_i_1_max__ = L;
            for (size_t i_0__ = 0; i_0__ < GL1_i_0_max__; ++i_0__) {
                for (size_t i_1__ = 0; i_1__ < GL1_i_1_max__; ++i_1__) {
                    check_greater_or_equal(function__, "GL1[i_0__][i_1__]", GL1[i_0__][i_1__], 0);
                    check_less_or_equal(function__, "GL1[i_0__][i_1__]", GL1[i_0__][i_1__], 1);
                }
            }
            current_statement_begin__ = 24;
            validate_non_negative_index("GL2", "N", N);
            validate_non_negative_index("GL2", "L", L);
            context__.validate_dims("data initialization", "GL2", "double", context__.to_vec(N,L));
            GL2 = std::vector<std::vector<double> >(N, std::vector<double>(L, double(0)));
            vals_r__ = context__.vals_r("GL2");
            pos__ = 0;
            size_t GL2_k_0_max__ = N;
            size_t GL2_k_1_max__ = L;
            for (size_t k_1__ = 0; k_1__ < GL2_k_1_max__; ++k_1__) {
                for (size_t k_0__ = 0; k_0__ < GL2_k_0_max__; ++k_0__) {
                    GL2[k_0__][k_1__] = vals_r__[pos__++];
                }
            }
            size_t GL2_i_0_max__ = N;
            size_t GL2_i_1_max__ = L;
            for (size_t i_0__ = 0; i_0__ < GL2_i_0_max__; ++i_0__) {
                for (size_t i_1__ = 0; i_1__ < GL2_i_1_max__; ++i_1__) {
                    check_greater_or_equal(function__, "GL2[i_0__][i_1__]", GL2[i_0__][i_1__], 0);
                    check_less_or_equal(function__, "GL2[i_0__][i_1__]", GL2[i_0__][i_1__], 1);
                }
            }
            current_statement_begin__ = 25;
            validate_non_negative_index("ploidy", "N", N);
            context__.validate_dims("data initialization", "ploidy", "double", context__.to_vec(N));
            ploidy = std::vector<double>(N, double(0));
            vals_r__ = context__.vals_r("ploidy");
            pos__ = 0;
            size_t ploidy_k_0_max__ = N;
            for (size_t k_0__ = 0; k_0__ < ploidy_k_0_max__; ++k_0__) {
                ploidy[k_0__] = vals_r__[pos__++];
            }
            size_t ploidy_i_0_max__ = N;
            for (size_t i_0__ = 0; i_0__ < ploidy_i_0_max__; ++i_0__) {
                check_greater_or_equal(function__, "ploidy[i_0__]", ploidy[i_0__], 0);
                check_less_or_equal(function__, "ploidy[i_0__]", ploidy[i_0__], 2);
            }
            current_statement_begin__ = 26;
            validate_non_negative_index("pids", "N", N);
            context__.validate_dims("data initialization", "pids", "int", context__.to_vec(N));
            pids = std::vector<int>(N, int(0));
            vals_i__ = context__.vals_i("pids");
            pos__ = 0;
            size_t pids_k_0_max__ = N;
            for (size_t k_0__ = 0; k_0__ < pids_k_0_max__; ++k_0__) {
                pids[k_0__] = vals_i__[pos__++];
            }
            size_t pids_i_0_max__ = N;
            for (size_t i_0__ = 0; i_0__ < pids_i_0_max__; ++i_0__) {
                check_greater_or_equal(function__, "pids[i_0__]", pids[i_0__], 0);
            }
            // initialize transformed data variables
            // execute transformed data statements
            // validate transformed data
            // validate, set parameter ranges
            num_params_r__ = 0U;
            param_ranges_i__.clear();
            current_statement_begin__ = 30;
            validate_non_negative_index("P", "J", J);
            validate_non_negative_index("P", "L", L);
            num_params_r__ += ((1 * J) * L);
        } catch (const std::exception& e) {
            stan::lang::rethrow_located(e, current_statement_begin__, prog_reader__());
            // Next line prevents compiler griping about no return
            throw std::runtime_error("*** IF YOU SEE THIS, PLEASE REPORT A BUG ***");
        }
    }
    ~model_popp_gl_mix() { }
    void transform_inits(const stan::io::var_context& context__,
                         std::vector<int>& params_i__,
                         std::vector<double>& params_r__,
                         std::ostream* pstream__) const {
        typedef double local_scalar_t__;
        stan::io::writer<double> writer__(params_r__, params_i__);
        size_t pos__;
        (void) pos__; // dummy call to supress warning
        std::vector<double> vals_r__;
        std::vector<int> vals_i__;
        current_statement_begin__ = 30;
        if (!(context__.contains_r("P")))
            stan::lang::rethrow_located(std::runtime_error(std::string("Variable P missing")), current_statement_begin__, prog_reader__());
        vals_r__ = context__.vals_r("P");
        pos__ = 0U;
        validate_non_negative_index("P", "J", J);
        validate_non_negative_index("P", "L", L);
        context__.validate_dims("parameter initialization", "P", "double", context__.to_vec(J,L));
        std::vector<std::vector<double> > P(J, std::vector<double>(L, double(0)));
        size_t P_k_0_max__ = J;
        size_t P_k_1_max__ = L;
        for (size_t k_1__ = 0; k_1__ < P_k_1_max__; ++k_1__) {
            for (size_t k_0__ = 0; k_0__ < P_k_0_max__; ++k_0__) {
                P[k_0__][k_1__] = vals_r__[pos__++];
            }
        }
        size_t P_i_0_max__ = J;
        size_t P_i_1_max__ = L;
        for (size_t i_0__ = 0; i_0__ < P_i_0_max__; ++i_0__) {
            for (size_t i_1__ = 0; i_1__ < P_i_1_max__; ++i_1__) {
                try {
                    writer__.scalar_lub_unconstrain(0, 1, P[i_0__][i_1__]);
                } catch (const std::exception& e) {
                    stan::lang::rethrow_located(std::runtime_error(std::string("Error transforming variable P: ") + e.what()), current_statement_begin__, prog_reader__());
                }
            }
        }
        params_r__ = writer__.data_r();
        params_i__ = writer__.data_i();
    }
    void transform_inits(const stan::io::var_context& context,
                         Eigen::Matrix<double, Eigen::Dynamic, 1>& params_r,
                         std::ostream* pstream__) const {
      std::vector<double> params_r_vec;
      std::vector<int> params_i_vec;
      transform_inits(context, params_i_vec, params_r_vec, pstream__);
      params_r.resize(params_r_vec.size());
      for (int i = 0; i < params_r.size(); ++i)
        params_r(i) = params_r_vec[i];
    }
    template <bool propto__, bool jacobian__, typename T__>
    T__ log_prob(std::vector<T__>& params_r__,
                 std::vector<int>& params_i__,
                 std::ostream* pstream__ = 0) const {
        typedef T__ local_scalar_t__;
        local_scalar_t__ DUMMY_VAR__(std::numeric_limits<double>::quiet_NaN());
        (void) DUMMY_VAR__;  // dummy to suppress unused var warning
        T__ lp__(0.0);
        stan::math::accumulator<T__> lp_accum__;
        try {
            stan::io::reader<local_scalar_t__> in__(params_r__, params_i__);
            // model parameters
            current_statement_begin__ = 30;
            std::vector<std::vector<local_scalar_t__> > P;
            size_t P_d_0_max__ = J;
            size_t P_d_1_max__ = L;
            P.resize(P_d_0_max__);
            for (size_t d_0__ = 0; d_0__ < P_d_0_max__; ++d_0__) {
                P[d_0__].reserve(P_d_1_max__);
                for (size_t d_1__ = 0; d_1__ < P_d_1_max__; ++d_1__) {
                    if (jacobian__)
                        P[d_0__].push_back(in__.scalar_lub_constrain(0, 1, lp__));
                    else
                        P[d_0__].push_back(in__.scalar_lub_constrain(0, 1));
                }
            }
            // model body
            current_statement_begin__ = 35;
            for (int i = 1; i <= L; ++i) {
                current_statement_begin__ = 36;
                for (int j = 1; j <= N; ++j) {
                    current_statement_begin__ = 38;
                    lp_accum__.add(calc_lik(get_base1(get_base1(GL0, j, "GL0", 1), i, "GL0", 2), get_base1(get_base1(GL1, j, "GL1", 1), i, "GL1", 2), get_base1(get_base1(GL2, j, "GL2", 1), i, "GL2", 2), get_base1(get_base1(P, get_base1(pids, j, "pids", 1), "P", 1), i, "P", 2), get_base1(ploidy, j, "ploidy", 1), pstream__));
                }
            }
            current_statement_begin__ = 41;
            for (int i = 1; i <= L; ++i) {
                current_statement_begin__ = 42;
                for (int j = 1; j <= J; ++j) {
                    current_statement_begin__ = 44;
                    lp_accum__.add(beta_log(get_base1(get_base1(P, j, "P", 1), i, "P", 2), 0.5, 0.5));
                }
            }
        } catch (const std::exception& e) {
            stan::lang::rethrow_located(e, current_statement_begin__, prog_reader__());
            // Next line prevents compiler griping about no return
            throw std::runtime_error("*** IF YOU SEE THIS, PLEASE REPORT A BUG ***");
        }
        lp_accum__.add(lp__);
        return lp_accum__.sum();
    } // log_prob()
    template <bool propto, bool jacobian, typename T_>
    T_ log_prob(Eigen::Matrix<T_,Eigen::Dynamic,1>& params_r,
               std::ostream* pstream = 0) const {
      std::vector<T_> vec_params_r;
      vec_params_r.reserve(params_r.size());
      for (int i = 0; i < params_r.size(); ++i)
        vec_params_r.push_back(params_r(i));
      std::vector<int> vec_params_i;
      return log_prob<propto,jacobian,T_>(vec_params_r, vec_params_i, pstream);
    }
    void get_param_names(std::vector<std::string>& names__) const {
        names__.resize(0);
        names__.push_back("P");
    }
    void get_dims(std::vector<std::vector<size_t> >& dimss__) const {
        dimss__.resize(0);
        std::vector<size_t> dims__;
        dims__.resize(0);
        dims__.push_back(J);
        dims__.push_back(L);
        dimss__.push_back(dims__);
    }
    template <typename RNG>
    void write_array(RNG& base_rng__,
                     std::vector<double>& params_r__,
                     std::vector<int>& params_i__,
                     std::vector<double>& vars__,
                     bool include_tparams__ = true,
                     bool include_gqs__ = true,
                     std::ostream* pstream__ = 0) const {
        typedef double local_scalar_t__;
        vars__.resize(0);
        stan::io::reader<local_scalar_t__> in__(params_r__, params_i__);
        static const char* function__ = "model_popp_gl_mix_namespace::write_array";
        (void) function__;  // dummy to suppress unused var warning
        // read-transform, write parameters
        std::vector<std::vector<double> > P;
        size_t P_d_0_max__ = J;
        size_t P_d_1_max__ = L;
        P.resize(P_d_0_max__);
        for (size_t d_0__ = 0; d_0__ < P_d_0_max__; ++d_0__) {
            P[d_0__].reserve(P_d_1_max__);
            for (size_t d_1__ = 0; d_1__ < P_d_1_max__; ++d_1__) {
                P[d_0__].push_back(in__.scalar_lub_constrain(0, 1));
            }
        }
        size_t P_k_0_max__ = J;
        size_t P_k_1_max__ = L;
        for (size_t k_1__ = 0; k_1__ < P_k_1_max__; ++k_1__) {
            for (size_t k_0__ = 0; k_0__ < P_k_0_max__; ++k_0__) {
                vars__.push_back(P[k_0__][k_1__]);
            }
        }
        double lp__ = 0.0;
        (void) lp__;  // dummy to suppress unused var warning
        stan::math::accumulator<double> lp_accum__;
        local_scalar_t__ DUMMY_VAR__(std::numeric_limits<double>::quiet_NaN());
        (void) DUMMY_VAR__;  // suppress unused var warning
        if (!include_tparams__ && !include_gqs__) return;
        try {
            if (!include_gqs__ && !include_tparams__) return;
            if (!include_gqs__) return;
        } catch (const std::exception& e) {
            stan::lang::rethrow_located(e, current_statement_begin__, prog_reader__());
            // Next line prevents compiler griping about no return
            throw std::runtime_error("*** IF YOU SEE THIS, PLEASE REPORT A BUG ***");
        }
    }
    template <typename RNG>
    void write_array(RNG& base_rng,
                     Eigen::Matrix<double,Eigen::Dynamic,1>& params_r,
                     Eigen::Matrix<double,Eigen::Dynamic,1>& vars,
                     bool include_tparams = true,
                     bool include_gqs = true,
                     std::ostream* pstream = 0) const {
      std::vector<double> params_r_vec(params_r.size());
      for (int i = 0; i < params_r.size(); ++i)
        params_r_vec[i] = params_r(i);
      std::vector<double> vars_vec;
      std::vector<int> params_i_vec;
      write_array(base_rng, params_r_vec, params_i_vec, vars_vec, include_tparams, include_gqs, pstream);
      vars.resize(vars_vec.size());
      for (int i = 0; i < vars.size(); ++i)
        vars(i) = vars_vec[i];
    }
    std::string model_name() const {
        return "model_popp_gl_mix";
    }
    void constrained_param_names(std::vector<std::string>& param_names__,
                                 bool include_tparams__ = true,
                                 bool include_gqs__ = true) const {
        std::stringstream param_name_stream__;
        size_t P_k_0_max__ = J;
        size_t P_k_1_max__ = L;
        for (size_t k_1__ = 0; k_1__ < P_k_1_max__; ++k_1__) {
            for (size_t k_0__ = 0; k_0__ < P_k_0_max__; ++k_0__) {
                param_name_stream__.str(std::string());
                param_name_stream__ << "P" << '.' << k_0__ + 1 << '.' << k_1__ + 1;
                param_names__.push_back(param_name_stream__.str());
            }
        }
        if (!include_gqs__ && !include_tparams__) return;
        if (include_tparams__) {
        }
        if (!include_gqs__) return;
    }
    void unconstrained_param_names(std::vector<std::string>& param_names__,
                                   bool include_tparams__ = true,
                                   bool include_gqs__ = true) const {
        std::stringstream param_name_stream__;
        size_t P_k_0_max__ = J;
        size_t P_k_1_max__ = L;
        for (size_t k_1__ = 0; k_1__ < P_k_1_max__; ++k_1__) {
            for (size_t k_0__ = 0; k_0__ < P_k_0_max__; ++k_0__) {
                param_name_stream__.str(std::string());
                param_name_stream__ << "P" << '.' << k_0__ + 1 << '.' << k_1__ + 1;
                param_names__.push_back(param_name_stream__.str());
            }
        }
        if (!include_gqs__ && !include_tparams__) return;
        if (include_tparams__) {
        }
        if (!include_gqs__) return;
    }
}; // model
}  // namespace
typedef model_popp_gl_mix_namespace::model_popp_gl_mix stan_model;
#ifndef USING_R
stan::model::model_base& new_model(
        stan::io::var_context& data_context,
        unsigned int seed,
        std::ostream* msg_stream) {
  stan_model* m = new stan_model(data_context, seed, msg_stream);
  return *m;
}
#endif
#endif
