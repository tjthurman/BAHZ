/*
    stantoolstest is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    stantoolstest is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with stantoolstest.  If not, see <http://www.gnu.org/licenses/>.
*/
#ifndef MODELS_HPP
#define MODELS_HPP
#define STAN__SERVICES__COMMAND_HPP
#include <rstan/rstaninc.hpp>
// Code generated by Stan version 2.17.0

#include <stan/model/model_header.hpp>

namespace model_binom_free_none_namespace {

using std::istream;
using std::string;
using std::stringstream;
using std::vector;
using stan::io::dump;
using stan::math::lgamma;
using stan::model::prob_grad;
using namespace stan::math;

typedef Eigen::Matrix<double,Eigen::Dynamic,1> vector_d;
typedef Eigen::Matrix<double,1,Eigen::Dynamic> row_vector_d;
typedef Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic> matrix_d;

static int current_statement_begin__;

stan::io::program_reader prog_reader__() {
    stan::io::program_reader reader;
    reader.add_event(0, 0, "start", "model_binom_free_none");
    reader.add_event(48, 48, "end", "model_binom_free_none");
    return reader;
}

#include <meta_header.hpp>
 class model_binom_free_none : public prob_grad {
private:
    int N;
    vector<int> nFocalAllele;
    vector<int> nTotalAlleles;
    vector<double> transectDist;
public:
    model_binom_free_none(stan::io::var_context& context__,
        std::ostream* pstream__ = 0)
        : prob_grad(0) {
        ctor_body(context__, 0, pstream__);
    }

    model_binom_free_none(stan::io::var_context& context__,
        unsigned int random_seed__,
        std::ostream* pstream__ = 0)
        : prob_grad(0) {
        ctor_body(context__, random_seed__, pstream__);
    }

    void ctor_body(stan::io::var_context& context__,
                   unsigned int random_seed__,
                   std::ostream* pstream__) {
        boost::ecuyer1988 base_rng__ =
          stan::services::util::create_rng(random_seed__, 0);
        (void) base_rng__;  // suppress unused var warning

        current_statement_begin__ = -1;

        static const char* function__ = "model_binom_free_none_namespace::model_binom_free_none";
        (void) function__;  // dummy to suppress unused var warning
        size_t pos__;
        (void) pos__;  // dummy to suppress unused var warning
        std::vector<int> vals_i__;
        std::vector<double> vals_r__;
        double DUMMY_VAR__(std::numeric_limits<double>::quiet_NaN());
        (void) DUMMY_VAR__;  // suppress unused var warning

        // initialize member variables
        try {
            current_statement_begin__ = 6;
            context__.validate_dims("data initialization", "N", "int", context__.to_vec());
            N = int(0);
            vals_i__ = context__.vals_i("N");
            pos__ = 0;
            N = vals_i__[pos__++];
            current_statement_begin__ = 7;
            validate_non_negative_index("nFocalAllele", "N", N);
            context__.validate_dims("data initialization", "nFocalAllele", "int", context__.to_vec(N));
            validate_non_negative_index("nFocalAllele", "N", N);
            nFocalAllele = std::vector<int>(N,int(0));
            vals_i__ = context__.vals_i("nFocalAllele");
            pos__ = 0;
            size_t nFocalAllele_limit_0__ = N;
            for (size_t i_0__ = 0; i_0__ < nFocalAllele_limit_0__; ++i_0__) {
                nFocalAllele[i_0__] = vals_i__[pos__++];
            }
            current_statement_begin__ = 8;
            validate_non_negative_index("nTotalAlleles", "N", N);
            context__.validate_dims("data initialization", "nTotalAlleles", "int", context__.to_vec(N));
            validate_non_negative_index("nTotalAlleles", "N", N);
            nTotalAlleles = std::vector<int>(N,int(0));
            vals_i__ = context__.vals_i("nTotalAlleles");
            pos__ = 0;
            size_t nTotalAlleles_limit_0__ = N;
            for (size_t i_0__ = 0; i_0__ < nTotalAlleles_limit_0__; ++i_0__) {
                nTotalAlleles[i_0__] = vals_i__[pos__++];
            }
            current_statement_begin__ = 9;
            validate_non_negative_index("transectDist", "N", N);
            context__.validate_dims("data initialization", "transectDist", "double", context__.to_vec(N));
            validate_non_negative_index("transectDist", "N", N);
            transectDist = std::vector<double>(N,double(0));
            vals_r__ = context__.vals_r("transectDist");
            pos__ = 0;
            size_t transectDist_limit_0__ = N;
            for (size_t i_0__ = 0; i_0__ < transectDist_limit_0__; ++i_0__) {
                transectDist[i_0__] = vals_r__[pos__++];
            }

            // validate, data variables
            current_statement_begin__ = 6;
            check_greater_or_equal(function__,"N",N,1);
            current_statement_begin__ = 7;
            current_statement_begin__ = 8;
            current_statement_begin__ = 9;
            // initialize data variables


            // validate transformed data

            // validate, set parameter ranges
            num_params_r__ = 0U;
            param_ranges_i__.clear();
            current_statement_begin__ = 12;
            ++num_params_r__;
            current_statement_begin__ = 13;
            ++num_params_r__;
            current_statement_begin__ = 14;
            ++num_params_r__;
            current_statement_begin__ = 15;
            ++num_params_r__;
        } catch (const std::exception& e) {
            stan::lang::rethrow_located(e, current_statement_begin__, prog_reader__());
            // Next line prevents compiler griping about no return
            throw std::runtime_error("*** IF YOU SEE THIS, PLEASE REPORT A BUG ***");
        }
    }

    ~model_binom_free_none() { }


    void transform_inits(const stan::io::var_context& context__,
                         std::vector<int>& params_i__,
                         std::vector<double>& params_r__,
                         std::ostream* pstream__) const {
        stan::io::writer<double> writer__(params_r__,params_i__);
        size_t pos__;
        (void) pos__; // dummy call to supress warning
        std::vector<double> vals_r__;
        std::vector<int> vals_i__;

        if (!(context__.contains_r("center")))
            throw std::runtime_error("variable center missing");
        vals_r__ = context__.vals_r("center");
        pos__ = 0U;
        context__.validate_dims("initialization", "center", "double", context__.to_vec());
        double center(0);
        center = vals_r__[pos__++];
        try {
            writer__.scalar_lb_unconstrain(0,center);
        } catch (const std::exception& e) { 
            throw std::runtime_error(std::string("Error transforming variable center: ") + e.what());
        }

        if (!(context__.contains_r("width")))
            throw std::runtime_error("variable width missing");
        vals_r__ = context__.vals_r("width");
        pos__ = 0U;
        context__.validate_dims("initialization", "width", "double", context__.to_vec());
        double width(0);
        width = vals_r__[pos__++];
        try {
            writer__.scalar_lb_unconstrain(0,width);
        } catch (const std::exception& e) { 
            throw std::runtime_error(std::string("Error transforming variable width: ") + e.what());
        }

        if (!(context__.contains_r("pmin")))
            throw std::runtime_error("variable pmin missing");
        vals_r__ = context__.vals_r("pmin");
        pos__ = 0U;
        context__.validate_dims("initialization", "pmin", "double", context__.to_vec());
        double pmin(0);
        pmin = vals_r__[pos__++];
        try {
            writer__.scalar_lub_unconstrain(0,1,pmin);
        } catch (const std::exception& e) { 
            throw std::runtime_error(std::string("Error transforming variable pmin: ") + e.what());
        }

        if (!(context__.contains_r("pmax")))
            throw std::runtime_error("variable pmax missing");
        vals_r__ = context__.vals_r("pmax");
        pos__ = 0U;
        context__.validate_dims("initialization", "pmax", "double", context__.to_vec());
        double pmax(0);
        pmax = vals_r__[pos__++];
        try {
            writer__.scalar_lub_unconstrain(0,1,pmax);
        } catch (const std::exception& e) { 
            throw std::runtime_error(std::string("Error transforming variable pmax: ") + e.what());
        }

        params_r__ = writer__.data_r();
        params_i__ = writer__.data_i();
    }

    void transform_inits(const stan::io::var_context& context,
                         Eigen::Matrix<double,Eigen::Dynamic,1>& params_r,
                         std::ostream* pstream__) const {
      std::vector<double> params_r_vec;
      std::vector<int> params_i_vec;
      transform_inits(context, params_i_vec, params_r_vec, pstream__);
      params_r.resize(params_r_vec.size());
      for (int i = 0; i < params_r.size(); ++i)
        params_r(i) = params_r_vec[i];
    }


    template <bool propto__, bool jacobian__, typename T__>
    T__ log_prob(vector<T__>& params_r__,
                 vector<int>& params_i__,
                 std::ostream* pstream__ = 0) const {

        T__ DUMMY_VAR__(std::numeric_limits<double>::quiet_NaN());
        (void) DUMMY_VAR__;  // suppress unused var warning

        T__ lp__(0.0);
        stan::math::accumulator<T__> lp_accum__;

        try {
            // model parameters
            stan::io::reader<T__> in__(params_r__,params_i__);

            T__ center;
            (void) center;  // dummy to suppress unused var warning
            if (jacobian__)
                center = in__.scalar_lb_constrain(0,lp__);
            else
                center = in__.scalar_lb_constrain(0);

            T__ width;
            (void) width;  // dummy to suppress unused var warning
            if (jacobian__)
                width = in__.scalar_lb_constrain(0,lp__);
            else
                width = in__.scalar_lb_constrain(0);

            T__ pmin;
            (void) pmin;  // dummy to suppress unused var warning
            if (jacobian__)
                pmin = in__.scalar_lub_constrain(0,1,lp__);
            else
                pmin = in__.scalar_lub_constrain(0,1);

            T__ pmax;
            (void) pmax;  // dummy to suppress unused var warning
            if (jacobian__)
                pmax = in__.scalar_lub_constrain(0,1,lp__);
            else
                pmax = in__.scalar_lub_constrain(0,1);


            // transformed parameters



            // validate transformed parameters

            const char* function__ = "validate transformed params";
            (void) function__;  // dummy to suppress unused var warning

            // model body
            {
            current_statement_begin__ = 19;
            validate_non_negative_index("p", "N", N);
            Eigen::Matrix<T__,Eigen::Dynamic,1>  p(static_cast<Eigen::VectorXd::Index>(N));
            (void) p;  // dummy to suppress unused var warning

            stan::math::initialize(p, DUMMY_VAR__);
            stan::math::fill(p,DUMMY_VAR__);


            current_statement_begin__ = 20;
            lp_accum__.add(uniform_log<propto__>(pmax, 0.80000000000000004, 1));
            current_statement_begin__ = 21;
            lp_accum__.add(uniform_log<propto__>(pmin, 0, 0.20000000000000001));
            current_statement_begin__ = 22;
            lp_accum__.add(normal_log<propto__>(width, 50, 100));
            current_statement_begin__ = 23;
            lp_accum__.add(normal_log<propto__>(center, 350, 150));
            current_statement_begin__ = 24;
            for (int i = 1; i <= N; ++i) {

                current_statement_begin__ = 26;
                stan::math::assign(get_base1_lhs(p,i,"p",1), (pmin + ((pmax - pmin) * (exp(((4 * (get_base1(transectDist,i,"transectDist",1) - center)) / width)) / (1 + exp(((4 * (get_base1(transectDist,i,"transectDist",1) - center)) / width)))))));
            }
            current_statement_begin__ = 30;
            lp_accum__.add(binomial_log<propto__>(nFocalAllele, nTotalAlleles, p));
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
        names__.push_back("center");
        names__.push_back("width");
        names__.push_back("pmin");
        names__.push_back("pmax");
        names__.push_back("p");
        names__.push_back("dev");
        names__.push_back("log_lik");
        names__.push_back("y_rep");
    }


    void get_dims(std::vector<std::vector<size_t> >& dimss__) const {
        dimss__.resize(0);
        std::vector<size_t> dims__;
        dims__.resize(0);
        dimss__.push_back(dims__);
        dims__.resize(0);
        dimss__.push_back(dims__);
        dims__.resize(0);
        dimss__.push_back(dims__);
        dims__.resize(0);
        dimss__.push_back(dims__);
        dims__.resize(0);
        dims__.push_back(N);
        dimss__.push_back(dims__);
        dims__.resize(0);
        dimss__.push_back(dims__);
        dims__.resize(0);
        dims__.push_back(N);
        dimss__.push_back(dims__);
        dims__.resize(0);
        dims__.push_back(N);
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
        vars__.resize(0);
        stan::io::reader<double> in__(params_r__,params_i__);
        static const char* function__ = "model_binom_free_none_namespace::write_array";
        (void) function__;  // dummy to suppress unused var warning
        // read-transform, write parameters
        double center = in__.scalar_lb_constrain(0);
        double width = in__.scalar_lb_constrain(0);
        double pmin = in__.scalar_lub_constrain(0,1);
        double pmax = in__.scalar_lub_constrain(0,1);
        vars__.push_back(center);
        vars__.push_back(width);
        vars__.push_back(pmin);
        vars__.push_back(pmax);

        if (!include_tparams__) return;
        // declare and define transformed parameters
        double lp__ = 0.0;
        (void) lp__;  // dummy to suppress unused var warning
        stan::math::accumulator<double> lp_accum__;

        double DUMMY_VAR__(std::numeric_limits<double>::quiet_NaN());
        (void) DUMMY_VAR__;  // suppress unused var warning

        try {



            // validate transformed parameters

            // write transformed parameters

            if (!include_gqs__) return;
            // declare and define generated quantities
            current_statement_begin__ = 34;
            validate_non_negative_index("p", "N", N);
            vector_d p(static_cast<Eigen::VectorXd::Index>(N));
            (void) p;  // dummy to suppress unused var warning

            stan::math::initialize(p, std::numeric_limits<double>::quiet_NaN());
            stan::math::fill(p,DUMMY_VAR__);
            current_statement_begin__ = 35;
            double dev(0.0);
            (void) dev;  // dummy to suppress unused var warning

            stan::math::initialize(dev, std::numeric_limits<double>::quiet_NaN());
            stan::math::fill(dev,DUMMY_VAR__);
            current_statement_begin__ = 36;
            validate_non_negative_index("log_lik", "N", N);
            vector_d log_lik(static_cast<Eigen::VectorXd::Index>(N));
            (void) log_lik;  // dummy to suppress unused var warning

            stan::math::initialize(log_lik, std::numeric_limits<double>::quiet_NaN());
            stan::math::fill(log_lik,DUMMY_VAR__);
            current_statement_begin__ = 37;
            validate_non_negative_index("y_rep", "N", N);
            vector_d y_rep(static_cast<Eigen::VectorXd::Index>(N));
            (void) y_rep;  // dummy to suppress unused var warning

            stan::math::initialize(y_rep, std::numeric_limits<double>::quiet_NaN());
            stan::math::fill(y_rep,DUMMY_VAR__);


            current_statement_begin__ = 38;
            for (int i = 1; i <= N; ++i) {

                current_statement_begin__ = 40;
                stan::math::assign(get_base1_lhs(p,i,"p",1), (pmin + ((pmax - pmin) * (exp(((4 * (get_base1(transectDist,i,"transectDist",1) - center)) / width)) / (1 + exp(((4 * (get_base1(transectDist,i,"transectDist",1) - center)) / width)))))));
                current_statement_begin__ = 42;
                stan::math::assign(get_base1_lhs(log_lik,i,"log_lik",1), binomial_log(get_base1(nFocalAllele,i,"nFocalAllele",1),get_base1(nTotalAlleles,i,"nTotalAlleles",1),get_base1(p,i,"p",1)));
                current_statement_begin__ = 44;
                stan::math::assign(get_base1_lhs(y_rep,i,"y_rep",1), binomial_rng(get_base1(nTotalAlleles,i,"nTotalAlleles",1),get_base1(p,i,"p",1), base_rng__));
            }
            current_statement_begin__ = 47;
            stan::math::assign(dev, (-(2) * binomial_log(nFocalAllele,nTotalAlleles,p)));

            // validate generated quantities
            current_statement_begin__ = 34;
            current_statement_begin__ = 35;
            current_statement_begin__ = 36;
            current_statement_begin__ = 37;

            // write generated quantities
            for (int k_0__ = 0; k_0__ < N; ++k_0__) {
            vars__.push_back(p[k_0__]);
            }
        vars__.push_back(dev);
            for (int k_0__ = 0; k_0__ < N; ++k_0__) {
            vars__.push_back(log_lik[k_0__]);
            }
            for (int k_0__ = 0; k_0__ < N; ++k_0__) {
            vars__.push_back(y_rep[k_0__]);
            }

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
      write_array(base_rng,params_r_vec,params_i_vec,vars_vec,include_tparams,include_gqs,pstream);
      vars.resize(vars_vec.size());
      for (int i = 0; i < vars.size(); ++i)
        vars(i) = vars_vec[i];
    }

    static std::string model_name() {
        return "model_binom_free_none";
    }


    void constrained_param_names(std::vector<std::string>& param_names__,
                                 bool include_tparams__ = true,
                                 bool include_gqs__ = true) const {
        std::stringstream param_name_stream__;
        param_name_stream__.str(std::string());
        param_name_stream__ << "center";
        param_names__.push_back(param_name_stream__.str());
        param_name_stream__.str(std::string());
        param_name_stream__ << "width";
        param_names__.push_back(param_name_stream__.str());
        param_name_stream__.str(std::string());
        param_name_stream__ << "pmin";
        param_names__.push_back(param_name_stream__.str());
        param_name_stream__.str(std::string());
        param_name_stream__ << "pmax";
        param_names__.push_back(param_name_stream__.str());

        if (!include_gqs__ && !include_tparams__) return;

        if (!include_gqs__) return;
        for (int k_0__ = 1; k_0__ <= N; ++k_0__) {
            param_name_stream__.str(std::string());
            param_name_stream__ << "p" << '.' << k_0__;
            param_names__.push_back(param_name_stream__.str());
        }
        param_name_stream__.str(std::string());
        param_name_stream__ << "dev";
        param_names__.push_back(param_name_stream__.str());
        for (int k_0__ = 1; k_0__ <= N; ++k_0__) {
            param_name_stream__.str(std::string());
            param_name_stream__ << "log_lik" << '.' << k_0__;
            param_names__.push_back(param_name_stream__.str());
        }
        for (int k_0__ = 1; k_0__ <= N; ++k_0__) {
            param_name_stream__.str(std::string());
            param_name_stream__ << "y_rep" << '.' << k_0__;
            param_names__.push_back(param_name_stream__.str());
        }
    }


    void unconstrained_param_names(std::vector<std::string>& param_names__,
                                   bool include_tparams__ = true,
                                   bool include_gqs__ = true) const {
        std::stringstream param_name_stream__;
        param_name_stream__.str(std::string());
        param_name_stream__ << "center";
        param_names__.push_back(param_name_stream__.str());
        param_name_stream__.str(std::string());
        param_name_stream__ << "width";
        param_names__.push_back(param_name_stream__.str());
        param_name_stream__.str(std::string());
        param_name_stream__ << "pmin";
        param_names__.push_back(param_name_stream__.str());
        param_name_stream__.str(std::string());
        param_name_stream__ << "pmax";
        param_names__.push_back(param_name_stream__.str());

        if (!include_gqs__ && !include_tparams__) return;

        if (!include_gqs__) return;
        for (int k_0__ = 1; k_0__ <= N; ++k_0__) {
            param_name_stream__.str(std::string());
            param_name_stream__ << "p" << '.' << k_0__;
            param_names__.push_back(param_name_stream__.str());
        }
        param_name_stream__.str(std::string());
        param_name_stream__ << "dev";
        param_names__.push_back(param_name_stream__.str());
        for (int k_0__ = 1; k_0__ <= N; ++k_0__) {
            param_name_stream__.str(std::string());
            param_name_stream__ << "log_lik" << '.' << k_0__;
            param_names__.push_back(param_name_stream__.str());
        }
        for (int k_0__ = 1; k_0__ <= N; ++k_0__) {
            param_name_stream__.str(std::string());
            param_name_stream__ << "y_rep" << '.' << k_0__;
            param_names__.push_back(param_name_stream__.str());
        }
    }

}; // model

}

typedef model_binom_free_none_namespace::model_binom_free_none stan_model;


#endif
