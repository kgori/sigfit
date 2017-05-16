#include <Rcpp.h>
using namespace Rcpp ;
#include "include/models.hpp"

RCPP_MODULE(stan_fit4sigfit_mod) {


    class_<rstan::stan_fit<model_sigfit_namespace::model_sigfit, boost::random::ecuyer1988> >("model_sigfit")

    .constructor<SEXP,SEXP>()


    .method("call_sampler", &rstan::stan_fit<model_sigfit_namespace::model_sigfit, boost::random::ecuyer1988> ::call_sampler)
    .method("param_names", &rstan::stan_fit<model_sigfit_namespace::model_sigfit, boost::random::ecuyer1988> ::param_names)
    .method("param_names_oi", &rstan::stan_fit<model_sigfit_namespace::model_sigfit, boost::random::ecuyer1988> ::param_names_oi)
    .method("param_fnames_oi", &rstan::stan_fit<model_sigfit_namespace::model_sigfit, boost::random::ecuyer1988> ::param_fnames_oi)
    .method("param_dims", &rstan::stan_fit<model_sigfit_namespace::model_sigfit, boost::random::ecuyer1988> ::param_dims)
    .method("param_dims_oi", &rstan::stan_fit<model_sigfit_namespace::model_sigfit, boost::random::ecuyer1988> ::param_dims_oi)
    .method("update_param_oi", &rstan::stan_fit<model_sigfit_namespace::model_sigfit, boost::random::ecuyer1988> ::update_param_oi)
    .method("param_oi_tidx", &rstan::stan_fit<model_sigfit_namespace::model_sigfit, boost::random::ecuyer1988> ::param_oi_tidx)
    .method("grad_log_prob", &rstan::stan_fit<model_sigfit_namespace::model_sigfit, boost::random::ecuyer1988> ::grad_log_prob)
    .method("log_prob", &rstan::stan_fit<model_sigfit_namespace::model_sigfit, boost::random::ecuyer1988> ::log_prob)
    .method("unconstrain_pars", &rstan::stan_fit<model_sigfit_namespace::model_sigfit, boost::random::ecuyer1988> ::unconstrain_pars)
    .method("constrain_pars", &rstan::stan_fit<model_sigfit_namespace::model_sigfit, boost::random::ecuyer1988> ::constrain_pars)
    .method("num_pars_unconstrained", &rstan::stan_fit<model_sigfit_namespace::model_sigfit, boost::random::ecuyer1988> ::num_pars_unconstrained)
    .method("unconstrained_param_names", &rstan::stan_fit<model_sigfit_namespace::model_sigfit, boost::random::ecuyer1988> ::unconstrained_param_names)
    .method("constrained_param_names", &rstan::stan_fit<model_sigfit_namespace::model_sigfit, boost::random::ecuyer1988> ::constrained_param_names)
    ;
}
#include <Rcpp.h>
using namespace Rcpp ;
#include "include/models.hpp"

RCPP_MODULE(stan_fit4sigfit_emu_mod) {


    class_<rstan::stan_fit<model_sigfit_emu_namespace::model_sigfit_emu, boost::random::ecuyer1988> >("model_sigfit_emu")

    .constructor<SEXP,SEXP>()


    .method("call_sampler", &rstan::stan_fit<model_sigfit_emu_namespace::model_sigfit_emu, boost::random::ecuyer1988> ::call_sampler)
    .method("param_names", &rstan::stan_fit<model_sigfit_emu_namespace::model_sigfit_emu, boost::random::ecuyer1988> ::param_names)
    .method("param_names_oi", &rstan::stan_fit<model_sigfit_emu_namespace::model_sigfit_emu, boost::random::ecuyer1988> ::param_names_oi)
    .method("param_fnames_oi", &rstan::stan_fit<model_sigfit_emu_namespace::model_sigfit_emu, boost::random::ecuyer1988> ::param_fnames_oi)
    .method("param_dims", &rstan::stan_fit<model_sigfit_emu_namespace::model_sigfit_emu, boost::random::ecuyer1988> ::param_dims)
    .method("param_dims_oi", &rstan::stan_fit<model_sigfit_emu_namespace::model_sigfit_emu, boost::random::ecuyer1988> ::param_dims_oi)
    .method("update_param_oi", &rstan::stan_fit<model_sigfit_emu_namespace::model_sigfit_emu, boost::random::ecuyer1988> ::update_param_oi)
    .method("param_oi_tidx", &rstan::stan_fit<model_sigfit_emu_namespace::model_sigfit_emu, boost::random::ecuyer1988> ::param_oi_tidx)
    .method("grad_log_prob", &rstan::stan_fit<model_sigfit_emu_namespace::model_sigfit_emu, boost::random::ecuyer1988> ::grad_log_prob)
    .method("log_prob", &rstan::stan_fit<model_sigfit_emu_namespace::model_sigfit_emu, boost::random::ecuyer1988> ::log_prob)
    .method("unconstrain_pars", &rstan::stan_fit<model_sigfit_emu_namespace::model_sigfit_emu, boost::random::ecuyer1988> ::unconstrain_pars)
    .method("constrain_pars", &rstan::stan_fit<model_sigfit_emu_namespace::model_sigfit_emu, boost::random::ecuyer1988> ::constrain_pars)
    .method("num_pars_unconstrained", &rstan::stan_fit<model_sigfit_emu_namespace::model_sigfit_emu, boost::random::ecuyer1988> ::num_pars_unconstrained)
    .method("unconstrained_param_names", &rstan::stan_fit<model_sigfit_emu_namespace::model_sigfit_emu, boost::random::ecuyer1988> ::unconstrained_param_names)
    .method("constrained_param_names", &rstan::stan_fit<model_sigfit_emu_namespace::model_sigfit_emu, boost::random::ecuyer1988> ::constrained_param_names)
    ;
}
#include <Rcpp.h>
using namespace Rcpp ;
#include "include/models.hpp"

RCPP_MODULE(stan_fit4sigfit_emu_fixed_mod) {


    class_<rstan::stan_fit<model_sigfit_emu_fixed_namespace::model_sigfit_emu_fixed, boost::random::ecuyer1988> >("model_sigfit_emu_fixed")

    .constructor<SEXP,SEXP>()


    .method("call_sampler", &rstan::stan_fit<model_sigfit_emu_fixed_namespace::model_sigfit_emu_fixed, boost::random::ecuyer1988> ::call_sampler)
    .method("param_names", &rstan::stan_fit<model_sigfit_emu_fixed_namespace::model_sigfit_emu_fixed, boost::random::ecuyer1988> ::param_names)
    .method("param_names_oi", &rstan::stan_fit<model_sigfit_emu_fixed_namespace::model_sigfit_emu_fixed, boost::random::ecuyer1988> ::param_names_oi)
    .method("param_fnames_oi", &rstan::stan_fit<model_sigfit_emu_fixed_namespace::model_sigfit_emu_fixed, boost::random::ecuyer1988> ::param_fnames_oi)
    .method("param_dims", &rstan::stan_fit<model_sigfit_emu_fixed_namespace::model_sigfit_emu_fixed, boost::random::ecuyer1988> ::param_dims)
    .method("param_dims_oi", &rstan::stan_fit<model_sigfit_emu_fixed_namespace::model_sigfit_emu_fixed, boost::random::ecuyer1988> ::param_dims_oi)
    .method("update_param_oi", &rstan::stan_fit<model_sigfit_emu_fixed_namespace::model_sigfit_emu_fixed, boost::random::ecuyer1988> ::update_param_oi)
    .method("param_oi_tidx", &rstan::stan_fit<model_sigfit_emu_fixed_namespace::model_sigfit_emu_fixed, boost::random::ecuyer1988> ::param_oi_tidx)
    .method("grad_log_prob", &rstan::stan_fit<model_sigfit_emu_fixed_namespace::model_sigfit_emu_fixed, boost::random::ecuyer1988> ::grad_log_prob)
    .method("log_prob", &rstan::stan_fit<model_sigfit_emu_fixed_namespace::model_sigfit_emu_fixed, boost::random::ecuyer1988> ::log_prob)
    .method("unconstrain_pars", &rstan::stan_fit<model_sigfit_emu_fixed_namespace::model_sigfit_emu_fixed, boost::random::ecuyer1988> ::unconstrain_pars)
    .method("constrain_pars", &rstan::stan_fit<model_sigfit_emu_fixed_namespace::model_sigfit_emu_fixed, boost::random::ecuyer1988> ::constrain_pars)
    .method("num_pars_unconstrained", &rstan::stan_fit<model_sigfit_emu_fixed_namespace::model_sigfit_emu_fixed, boost::random::ecuyer1988> ::num_pars_unconstrained)
    .method("unconstrained_param_names", &rstan::stan_fit<model_sigfit_emu_fixed_namespace::model_sigfit_emu_fixed, boost::random::ecuyer1988> ::unconstrained_param_names)
    .method("constrained_param_names", &rstan::stan_fit<model_sigfit_emu_fixed_namespace::model_sigfit_emu_fixed, boost::random::ecuyer1988> ::constrained_param_names)
    ;
}
#include <Rcpp.h>
using namespace Rcpp ;
#include "include/models.hpp"

RCPP_MODULE(stan_fit4sigfit_hier_mod) {


    class_<rstan::stan_fit<model_sigfit_hier_namespace::model_sigfit_hier, boost::random::ecuyer1988> >("model_sigfit_hier")

    .constructor<SEXP,SEXP>()


    .method("call_sampler", &rstan::stan_fit<model_sigfit_hier_namespace::model_sigfit_hier, boost::random::ecuyer1988> ::call_sampler)
    .method("param_names", &rstan::stan_fit<model_sigfit_hier_namespace::model_sigfit_hier, boost::random::ecuyer1988> ::param_names)
    .method("param_names_oi", &rstan::stan_fit<model_sigfit_hier_namespace::model_sigfit_hier, boost::random::ecuyer1988> ::param_names_oi)
    .method("param_fnames_oi", &rstan::stan_fit<model_sigfit_hier_namespace::model_sigfit_hier, boost::random::ecuyer1988> ::param_fnames_oi)
    .method("param_dims", &rstan::stan_fit<model_sigfit_hier_namespace::model_sigfit_hier, boost::random::ecuyer1988> ::param_dims)
    .method("param_dims_oi", &rstan::stan_fit<model_sigfit_hier_namespace::model_sigfit_hier, boost::random::ecuyer1988> ::param_dims_oi)
    .method("update_param_oi", &rstan::stan_fit<model_sigfit_hier_namespace::model_sigfit_hier, boost::random::ecuyer1988> ::update_param_oi)
    .method("param_oi_tidx", &rstan::stan_fit<model_sigfit_hier_namespace::model_sigfit_hier, boost::random::ecuyer1988> ::param_oi_tidx)
    .method("grad_log_prob", &rstan::stan_fit<model_sigfit_hier_namespace::model_sigfit_hier, boost::random::ecuyer1988> ::grad_log_prob)
    .method("log_prob", &rstan::stan_fit<model_sigfit_hier_namespace::model_sigfit_hier, boost::random::ecuyer1988> ::log_prob)
    .method("unconstrain_pars", &rstan::stan_fit<model_sigfit_hier_namespace::model_sigfit_hier, boost::random::ecuyer1988> ::unconstrain_pars)
    .method("constrain_pars", &rstan::stan_fit<model_sigfit_hier_namespace::model_sigfit_hier, boost::random::ecuyer1988> ::constrain_pars)
    .method("num_pars_unconstrained", &rstan::stan_fit<model_sigfit_hier_namespace::model_sigfit_hier, boost::random::ecuyer1988> ::num_pars_unconstrained)
    .method("unconstrained_param_names", &rstan::stan_fit<model_sigfit_hier_namespace::model_sigfit_hier, boost::random::ecuyer1988> ::unconstrained_param_names)
    .method("constrained_param_names", &rstan::stan_fit<model_sigfit_hier_namespace::model_sigfit_hier, boost::random::ecuyer1988> ::constrained_param_names)
    ;
}
#include <Rcpp.h>
using namespace Rcpp ;
#include "include/models.hpp"

RCPP_MODULE(stan_fit4sigfit_multi_mod) {


    class_<rstan::stan_fit<model_sigfit_multi_namespace::model_sigfit_multi, boost::random::ecuyer1988> >("model_sigfit_multi")

    .constructor<SEXP,SEXP>()


    .method("call_sampler", &rstan::stan_fit<model_sigfit_multi_namespace::model_sigfit_multi, boost::random::ecuyer1988> ::call_sampler)
    .method("param_names", &rstan::stan_fit<model_sigfit_multi_namespace::model_sigfit_multi, boost::random::ecuyer1988> ::param_names)
    .method("param_names_oi", &rstan::stan_fit<model_sigfit_multi_namespace::model_sigfit_multi, boost::random::ecuyer1988> ::param_names_oi)
    .method("param_fnames_oi", &rstan::stan_fit<model_sigfit_multi_namespace::model_sigfit_multi, boost::random::ecuyer1988> ::param_fnames_oi)
    .method("param_dims", &rstan::stan_fit<model_sigfit_multi_namespace::model_sigfit_multi, boost::random::ecuyer1988> ::param_dims)
    .method("param_dims_oi", &rstan::stan_fit<model_sigfit_multi_namespace::model_sigfit_multi, boost::random::ecuyer1988> ::param_dims_oi)
    .method("update_param_oi", &rstan::stan_fit<model_sigfit_multi_namespace::model_sigfit_multi, boost::random::ecuyer1988> ::update_param_oi)
    .method("param_oi_tidx", &rstan::stan_fit<model_sigfit_multi_namespace::model_sigfit_multi, boost::random::ecuyer1988> ::param_oi_tidx)
    .method("grad_log_prob", &rstan::stan_fit<model_sigfit_multi_namespace::model_sigfit_multi, boost::random::ecuyer1988> ::grad_log_prob)
    .method("log_prob", &rstan::stan_fit<model_sigfit_multi_namespace::model_sigfit_multi, boost::random::ecuyer1988> ::log_prob)
    .method("unconstrain_pars", &rstan::stan_fit<model_sigfit_multi_namespace::model_sigfit_multi, boost::random::ecuyer1988> ::unconstrain_pars)
    .method("constrain_pars", &rstan::stan_fit<model_sigfit_multi_namespace::model_sigfit_multi, boost::random::ecuyer1988> ::constrain_pars)
    .method("num_pars_unconstrained", &rstan::stan_fit<model_sigfit_multi_namespace::model_sigfit_multi, boost::random::ecuyer1988> ::num_pars_unconstrained)
    .method("unconstrained_param_names", &rstan::stan_fit<model_sigfit_multi_namespace::model_sigfit_multi, boost::random::ecuyer1988> ::unconstrained_param_names)
    .method("constrained_param_names", &rstan::stan_fit<model_sigfit_multi_namespace::model_sigfit_multi, boost::random::ecuyer1988> ::constrained_param_names)
    ;
}
#include <Rcpp.h>
using namespace Rcpp ;
#include "include/models.hpp"

RCPP_MODULE(stan_fit4sigfit_nmf_mod) {


    class_<rstan::stan_fit<model_sigfit_nmf_namespace::model_sigfit_nmf, boost::random::ecuyer1988> >("model_sigfit_nmf")

    .constructor<SEXP,SEXP>()


    .method("call_sampler", &rstan::stan_fit<model_sigfit_nmf_namespace::model_sigfit_nmf, boost::random::ecuyer1988> ::call_sampler)
    .method("param_names", &rstan::stan_fit<model_sigfit_nmf_namespace::model_sigfit_nmf, boost::random::ecuyer1988> ::param_names)
    .method("param_names_oi", &rstan::stan_fit<model_sigfit_nmf_namespace::model_sigfit_nmf, boost::random::ecuyer1988> ::param_names_oi)
    .method("param_fnames_oi", &rstan::stan_fit<model_sigfit_nmf_namespace::model_sigfit_nmf, boost::random::ecuyer1988> ::param_fnames_oi)
    .method("param_dims", &rstan::stan_fit<model_sigfit_nmf_namespace::model_sigfit_nmf, boost::random::ecuyer1988> ::param_dims)
    .method("param_dims_oi", &rstan::stan_fit<model_sigfit_nmf_namespace::model_sigfit_nmf, boost::random::ecuyer1988> ::param_dims_oi)
    .method("update_param_oi", &rstan::stan_fit<model_sigfit_nmf_namespace::model_sigfit_nmf, boost::random::ecuyer1988> ::update_param_oi)
    .method("param_oi_tidx", &rstan::stan_fit<model_sigfit_nmf_namespace::model_sigfit_nmf, boost::random::ecuyer1988> ::param_oi_tidx)
    .method("grad_log_prob", &rstan::stan_fit<model_sigfit_nmf_namespace::model_sigfit_nmf, boost::random::ecuyer1988> ::grad_log_prob)
    .method("log_prob", &rstan::stan_fit<model_sigfit_nmf_namespace::model_sigfit_nmf, boost::random::ecuyer1988> ::log_prob)
    .method("unconstrain_pars", &rstan::stan_fit<model_sigfit_nmf_namespace::model_sigfit_nmf, boost::random::ecuyer1988> ::unconstrain_pars)
    .method("constrain_pars", &rstan::stan_fit<model_sigfit_nmf_namespace::model_sigfit_nmf, boost::random::ecuyer1988> ::constrain_pars)
    .method("num_pars_unconstrained", &rstan::stan_fit<model_sigfit_nmf_namespace::model_sigfit_nmf, boost::random::ecuyer1988> ::num_pars_unconstrained)
    .method("unconstrained_param_names", &rstan::stan_fit<model_sigfit_nmf_namespace::model_sigfit_nmf, boost::random::ecuyer1988> ::unconstrained_param_names)
    .method("constrained_param_names", &rstan::stan_fit<model_sigfit_nmf_namespace::model_sigfit_nmf, boost::random::ecuyer1988> ::constrained_param_names)
    ;
}
