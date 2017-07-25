#include <Rcpp.h>
using namespace Rcpp ;
#include "include/models.hpp"

RCPP_MODULE(stan_fit4sigfit_ext_emu_mod) {


    class_<rstan::stan_fit<model_sigfit_ext_emu_namespace::model_sigfit_ext_emu, boost::random::ecuyer1988> >("model_sigfit_ext_emu")

    .constructor<SEXP,SEXP>()


    .method("call_sampler", &rstan::stan_fit<model_sigfit_ext_emu_namespace::model_sigfit_ext_emu, boost::random::ecuyer1988> ::call_sampler)
    .method("param_names", &rstan::stan_fit<model_sigfit_ext_emu_namespace::model_sigfit_ext_emu, boost::random::ecuyer1988> ::param_names)
    .method("param_names_oi", &rstan::stan_fit<model_sigfit_ext_emu_namespace::model_sigfit_ext_emu, boost::random::ecuyer1988> ::param_names_oi)
    .method("param_fnames_oi", &rstan::stan_fit<model_sigfit_ext_emu_namespace::model_sigfit_ext_emu, boost::random::ecuyer1988> ::param_fnames_oi)
    .method("param_dims", &rstan::stan_fit<model_sigfit_ext_emu_namespace::model_sigfit_ext_emu, boost::random::ecuyer1988> ::param_dims)
    .method("param_dims_oi", &rstan::stan_fit<model_sigfit_ext_emu_namespace::model_sigfit_ext_emu, boost::random::ecuyer1988> ::param_dims_oi)
    .method("update_param_oi", &rstan::stan_fit<model_sigfit_ext_emu_namespace::model_sigfit_ext_emu, boost::random::ecuyer1988> ::update_param_oi)
    .method("param_oi_tidx", &rstan::stan_fit<model_sigfit_ext_emu_namespace::model_sigfit_ext_emu, boost::random::ecuyer1988> ::param_oi_tidx)
    .method("grad_log_prob", &rstan::stan_fit<model_sigfit_ext_emu_namespace::model_sigfit_ext_emu, boost::random::ecuyer1988> ::grad_log_prob)
    .method("log_prob", &rstan::stan_fit<model_sigfit_ext_emu_namespace::model_sigfit_ext_emu, boost::random::ecuyer1988> ::log_prob)
    .method("unconstrain_pars", &rstan::stan_fit<model_sigfit_ext_emu_namespace::model_sigfit_ext_emu, boost::random::ecuyer1988> ::unconstrain_pars)
    .method("constrain_pars", &rstan::stan_fit<model_sigfit_ext_emu_namespace::model_sigfit_ext_emu, boost::random::ecuyer1988> ::constrain_pars)
    .method("num_pars_unconstrained", &rstan::stan_fit<model_sigfit_ext_emu_namespace::model_sigfit_ext_emu, boost::random::ecuyer1988> ::num_pars_unconstrained)
    .method("unconstrained_param_names", &rstan::stan_fit<model_sigfit_ext_emu_namespace::model_sigfit_ext_emu, boost::random::ecuyer1988> ::unconstrained_param_names)
    .method("constrained_param_names", &rstan::stan_fit<model_sigfit_ext_emu_namespace::model_sigfit_ext_emu, boost::random::ecuyer1988> ::constrained_param_names)
    ;
}
#include <Rcpp.h>
using namespace Rcpp ;
#include "include/models.hpp"

RCPP_MODULE(stan_fit4sigfit_ext_nmf_mod) {


    class_<rstan::stan_fit<model_sigfit_ext_nmf_namespace::model_sigfit_ext_nmf, boost::random::ecuyer1988> >("model_sigfit_ext_nmf")

    .constructor<SEXP,SEXP>()


    .method("call_sampler", &rstan::stan_fit<model_sigfit_ext_nmf_namespace::model_sigfit_ext_nmf, boost::random::ecuyer1988> ::call_sampler)
    .method("param_names", &rstan::stan_fit<model_sigfit_ext_nmf_namespace::model_sigfit_ext_nmf, boost::random::ecuyer1988> ::param_names)
    .method("param_names_oi", &rstan::stan_fit<model_sigfit_ext_nmf_namespace::model_sigfit_ext_nmf, boost::random::ecuyer1988> ::param_names_oi)
    .method("param_fnames_oi", &rstan::stan_fit<model_sigfit_ext_nmf_namespace::model_sigfit_ext_nmf, boost::random::ecuyer1988> ::param_fnames_oi)
    .method("param_dims", &rstan::stan_fit<model_sigfit_ext_nmf_namespace::model_sigfit_ext_nmf, boost::random::ecuyer1988> ::param_dims)
    .method("param_dims_oi", &rstan::stan_fit<model_sigfit_ext_nmf_namespace::model_sigfit_ext_nmf, boost::random::ecuyer1988> ::param_dims_oi)
    .method("update_param_oi", &rstan::stan_fit<model_sigfit_ext_nmf_namespace::model_sigfit_ext_nmf, boost::random::ecuyer1988> ::update_param_oi)
    .method("param_oi_tidx", &rstan::stan_fit<model_sigfit_ext_nmf_namespace::model_sigfit_ext_nmf, boost::random::ecuyer1988> ::param_oi_tidx)
    .method("grad_log_prob", &rstan::stan_fit<model_sigfit_ext_nmf_namespace::model_sigfit_ext_nmf, boost::random::ecuyer1988> ::grad_log_prob)
    .method("log_prob", &rstan::stan_fit<model_sigfit_ext_nmf_namespace::model_sigfit_ext_nmf, boost::random::ecuyer1988> ::log_prob)
    .method("unconstrain_pars", &rstan::stan_fit<model_sigfit_ext_nmf_namespace::model_sigfit_ext_nmf, boost::random::ecuyer1988> ::unconstrain_pars)
    .method("constrain_pars", &rstan::stan_fit<model_sigfit_ext_nmf_namespace::model_sigfit_ext_nmf, boost::random::ecuyer1988> ::constrain_pars)
    .method("num_pars_unconstrained", &rstan::stan_fit<model_sigfit_ext_nmf_namespace::model_sigfit_ext_nmf, boost::random::ecuyer1988> ::num_pars_unconstrained)
    .method("unconstrained_param_names", &rstan::stan_fit<model_sigfit_ext_nmf_namespace::model_sigfit_ext_nmf, boost::random::ecuyer1988> ::unconstrained_param_names)
    .method("constrained_param_names", &rstan::stan_fit<model_sigfit_ext_nmf_namespace::model_sigfit_ext_nmf, boost::random::ecuyer1988> ::constrained_param_names)
    ;
}
#include <Rcpp.h>
using namespace Rcpp ;
#include "include/models.hpp"

RCPP_MODULE(stan_fit4sigfit_fit_emu_mod) {


    class_<rstan::stan_fit<model_sigfit_fit_emu_namespace::model_sigfit_fit_emu, boost::random::ecuyer1988> >("model_sigfit_fit_emu")

    .constructor<SEXP,SEXP>()


    .method("call_sampler", &rstan::stan_fit<model_sigfit_fit_emu_namespace::model_sigfit_fit_emu, boost::random::ecuyer1988> ::call_sampler)
    .method("param_names", &rstan::stan_fit<model_sigfit_fit_emu_namespace::model_sigfit_fit_emu, boost::random::ecuyer1988> ::param_names)
    .method("param_names_oi", &rstan::stan_fit<model_sigfit_fit_emu_namespace::model_sigfit_fit_emu, boost::random::ecuyer1988> ::param_names_oi)
    .method("param_fnames_oi", &rstan::stan_fit<model_sigfit_fit_emu_namespace::model_sigfit_fit_emu, boost::random::ecuyer1988> ::param_fnames_oi)
    .method("param_dims", &rstan::stan_fit<model_sigfit_fit_emu_namespace::model_sigfit_fit_emu, boost::random::ecuyer1988> ::param_dims)
    .method("param_dims_oi", &rstan::stan_fit<model_sigfit_fit_emu_namespace::model_sigfit_fit_emu, boost::random::ecuyer1988> ::param_dims_oi)
    .method("update_param_oi", &rstan::stan_fit<model_sigfit_fit_emu_namespace::model_sigfit_fit_emu, boost::random::ecuyer1988> ::update_param_oi)
    .method("param_oi_tidx", &rstan::stan_fit<model_sigfit_fit_emu_namespace::model_sigfit_fit_emu, boost::random::ecuyer1988> ::param_oi_tidx)
    .method("grad_log_prob", &rstan::stan_fit<model_sigfit_fit_emu_namespace::model_sigfit_fit_emu, boost::random::ecuyer1988> ::grad_log_prob)
    .method("log_prob", &rstan::stan_fit<model_sigfit_fit_emu_namespace::model_sigfit_fit_emu, boost::random::ecuyer1988> ::log_prob)
    .method("unconstrain_pars", &rstan::stan_fit<model_sigfit_fit_emu_namespace::model_sigfit_fit_emu, boost::random::ecuyer1988> ::unconstrain_pars)
    .method("constrain_pars", &rstan::stan_fit<model_sigfit_fit_emu_namespace::model_sigfit_fit_emu, boost::random::ecuyer1988> ::constrain_pars)
    .method("num_pars_unconstrained", &rstan::stan_fit<model_sigfit_fit_emu_namespace::model_sigfit_fit_emu, boost::random::ecuyer1988> ::num_pars_unconstrained)
    .method("unconstrained_param_names", &rstan::stan_fit<model_sigfit_fit_emu_namespace::model_sigfit_fit_emu, boost::random::ecuyer1988> ::unconstrained_param_names)
    .method("constrained_param_names", &rstan::stan_fit<model_sigfit_fit_emu_namespace::model_sigfit_fit_emu, boost::random::ecuyer1988> ::constrained_param_names)
    ;
}
#include <Rcpp.h>
using namespace Rcpp ;
#include "include/models.hpp"

RCPP_MODULE(stan_fit4sigfit_fit_nmf_mod) {


    class_<rstan::stan_fit<model_sigfit_fit_nmf_namespace::model_sigfit_fit_nmf, boost::random::ecuyer1988> >("model_sigfit_fit_nmf")

    .constructor<SEXP,SEXP>()


    .method("call_sampler", &rstan::stan_fit<model_sigfit_fit_nmf_namespace::model_sigfit_fit_nmf, boost::random::ecuyer1988> ::call_sampler)
    .method("param_names", &rstan::stan_fit<model_sigfit_fit_nmf_namespace::model_sigfit_fit_nmf, boost::random::ecuyer1988> ::param_names)
    .method("param_names_oi", &rstan::stan_fit<model_sigfit_fit_nmf_namespace::model_sigfit_fit_nmf, boost::random::ecuyer1988> ::param_names_oi)
    .method("param_fnames_oi", &rstan::stan_fit<model_sigfit_fit_nmf_namespace::model_sigfit_fit_nmf, boost::random::ecuyer1988> ::param_fnames_oi)
    .method("param_dims", &rstan::stan_fit<model_sigfit_fit_nmf_namespace::model_sigfit_fit_nmf, boost::random::ecuyer1988> ::param_dims)
    .method("param_dims_oi", &rstan::stan_fit<model_sigfit_fit_nmf_namespace::model_sigfit_fit_nmf, boost::random::ecuyer1988> ::param_dims_oi)
    .method("update_param_oi", &rstan::stan_fit<model_sigfit_fit_nmf_namespace::model_sigfit_fit_nmf, boost::random::ecuyer1988> ::update_param_oi)
    .method("param_oi_tidx", &rstan::stan_fit<model_sigfit_fit_nmf_namespace::model_sigfit_fit_nmf, boost::random::ecuyer1988> ::param_oi_tidx)
    .method("grad_log_prob", &rstan::stan_fit<model_sigfit_fit_nmf_namespace::model_sigfit_fit_nmf, boost::random::ecuyer1988> ::grad_log_prob)
    .method("log_prob", &rstan::stan_fit<model_sigfit_fit_nmf_namespace::model_sigfit_fit_nmf, boost::random::ecuyer1988> ::log_prob)
    .method("unconstrain_pars", &rstan::stan_fit<model_sigfit_fit_nmf_namespace::model_sigfit_fit_nmf, boost::random::ecuyer1988> ::unconstrain_pars)
    .method("constrain_pars", &rstan::stan_fit<model_sigfit_fit_nmf_namespace::model_sigfit_fit_nmf, boost::random::ecuyer1988> ::constrain_pars)
    .method("num_pars_unconstrained", &rstan::stan_fit<model_sigfit_fit_nmf_namespace::model_sigfit_fit_nmf, boost::random::ecuyer1988> ::num_pars_unconstrained)
    .method("unconstrained_param_names", &rstan::stan_fit<model_sigfit_fit_nmf_namespace::model_sigfit_fit_nmf, boost::random::ecuyer1988> ::unconstrained_param_names)
    .method("constrained_param_names", &rstan::stan_fit<model_sigfit_fit_nmf_namespace::model_sigfit_fit_nmf, boost::random::ecuyer1988> ::constrained_param_names)
    ;
}
#include <Rcpp.h>
using namespace Rcpp ;
#include "include/models.hpp"

RCPP_MODULE(stan_fit4sigfit_fit_nmf_hier_mod) {


    class_<rstan::stan_fit<model_sigfit_fit_nmf_hier_namespace::model_sigfit_fit_nmf_hier, boost::random::ecuyer1988> >("model_sigfit_fit_nmf_hier")

    .constructor<SEXP,SEXP>()


    .method("call_sampler", &rstan::stan_fit<model_sigfit_fit_nmf_hier_namespace::model_sigfit_fit_nmf_hier, boost::random::ecuyer1988> ::call_sampler)
    .method("param_names", &rstan::stan_fit<model_sigfit_fit_nmf_hier_namespace::model_sigfit_fit_nmf_hier, boost::random::ecuyer1988> ::param_names)
    .method("param_names_oi", &rstan::stan_fit<model_sigfit_fit_nmf_hier_namespace::model_sigfit_fit_nmf_hier, boost::random::ecuyer1988> ::param_names_oi)
    .method("param_fnames_oi", &rstan::stan_fit<model_sigfit_fit_nmf_hier_namespace::model_sigfit_fit_nmf_hier, boost::random::ecuyer1988> ::param_fnames_oi)
    .method("param_dims", &rstan::stan_fit<model_sigfit_fit_nmf_hier_namespace::model_sigfit_fit_nmf_hier, boost::random::ecuyer1988> ::param_dims)
    .method("param_dims_oi", &rstan::stan_fit<model_sigfit_fit_nmf_hier_namespace::model_sigfit_fit_nmf_hier, boost::random::ecuyer1988> ::param_dims_oi)
    .method("update_param_oi", &rstan::stan_fit<model_sigfit_fit_nmf_hier_namespace::model_sigfit_fit_nmf_hier, boost::random::ecuyer1988> ::update_param_oi)
    .method("param_oi_tidx", &rstan::stan_fit<model_sigfit_fit_nmf_hier_namespace::model_sigfit_fit_nmf_hier, boost::random::ecuyer1988> ::param_oi_tidx)
    .method("grad_log_prob", &rstan::stan_fit<model_sigfit_fit_nmf_hier_namespace::model_sigfit_fit_nmf_hier, boost::random::ecuyer1988> ::grad_log_prob)
    .method("log_prob", &rstan::stan_fit<model_sigfit_fit_nmf_hier_namespace::model_sigfit_fit_nmf_hier, boost::random::ecuyer1988> ::log_prob)
    .method("unconstrain_pars", &rstan::stan_fit<model_sigfit_fit_nmf_hier_namespace::model_sigfit_fit_nmf_hier, boost::random::ecuyer1988> ::unconstrain_pars)
    .method("constrain_pars", &rstan::stan_fit<model_sigfit_fit_nmf_hier_namespace::model_sigfit_fit_nmf_hier, boost::random::ecuyer1988> ::constrain_pars)
    .method("num_pars_unconstrained", &rstan::stan_fit<model_sigfit_fit_nmf_hier_namespace::model_sigfit_fit_nmf_hier, boost::random::ecuyer1988> ::num_pars_unconstrained)
    .method("unconstrained_param_names", &rstan::stan_fit<model_sigfit_fit_nmf_hier_namespace::model_sigfit_fit_nmf_hier, boost::random::ecuyer1988> ::unconstrained_param_names)
    .method("constrained_param_names", &rstan::stan_fit<model_sigfit_fit_nmf_hier_namespace::model_sigfit_fit_nmf_hier, boost::random::ecuyer1988> ::constrained_param_names)
    ;
}
#include <Rcpp.h>
using namespace Rcpp ;
#include "include/models.hpp"

RCPP_MODULE(stan_fit4sigfit_fitex_nmf_mod) {


    class_<rstan::stan_fit<model_sigfit_fitex_nmf_namespace::model_sigfit_fitex_nmf, boost::random::ecuyer1988> >("model_sigfit_fitex_nmf")

    .constructor<SEXP,SEXP>()


    .method("call_sampler", &rstan::stan_fit<model_sigfit_fitex_nmf_namespace::model_sigfit_fitex_nmf, boost::random::ecuyer1988> ::call_sampler)
    .method("param_names", &rstan::stan_fit<model_sigfit_fitex_nmf_namespace::model_sigfit_fitex_nmf, boost::random::ecuyer1988> ::param_names)
    .method("param_names_oi", &rstan::stan_fit<model_sigfit_fitex_nmf_namespace::model_sigfit_fitex_nmf, boost::random::ecuyer1988> ::param_names_oi)
    .method("param_fnames_oi", &rstan::stan_fit<model_sigfit_fitex_nmf_namespace::model_sigfit_fitex_nmf, boost::random::ecuyer1988> ::param_fnames_oi)
    .method("param_dims", &rstan::stan_fit<model_sigfit_fitex_nmf_namespace::model_sigfit_fitex_nmf, boost::random::ecuyer1988> ::param_dims)
    .method("param_dims_oi", &rstan::stan_fit<model_sigfit_fitex_nmf_namespace::model_sigfit_fitex_nmf, boost::random::ecuyer1988> ::param_dims_oi)
    .method("update_param_oi", &rstan::stan_fit<model_sigfit_fitex_nmf_namespace::model_sigfit_fitex_nmf, boost::random::ecuyer1988> ::update_param_oi)
    .method("param_oi_tidx", &rstan::stan_fit<model_sigfit_fitex_nmf_namespace::model_sigfit_fitex_nmf, boost::random::ecuyer1988> ::param_oi_tidx)
    .method("grad_log_prob", &rstan::stan_fit<model_sigfit_fitex_nmf_namespace::model_sigfit_fitex_nmf, boost::random::ecuyer1988> ::grad_log_prob)
    .method("log_prob", &rstan::stan_fit<model_sigfit_fitex_nmf_namespace::model_sigfit_fitex_nmf, boost::random::ecuyer1988> ::log_prob)
    .method("unconstrain_pars", &rstan::stan_fit<model_sigfit_fitex_nmf_namespace::model_sigfit_fitex_nmf, boost::random::ecuyer1988> ::unconstrain_pars)
    .method("constrain_pars", &rstan::stan_fit<model_sigfit_fitex_nmf_namespace::model_sigfit_fitex_nmf, boost::random::ecuyer1988> ::constrain_pars)
    .method("num_pars_unconstrained", &rstan::stan_fit<model_sigfit_fitex_nmf_namespace::model_sigfit_fitex_nmf, boost::random::ecuyer1988> ::num_pars_unconstrained)
    .method("unconstrained_param_names", &rstan::stan_fit<model_sigfit_fitex_nmf_namespace::model_sigfit_fitex_nmf, boost::random::ecuyer1988> ::unconstrained_param_names)
    .method("constrained_param_names", &rstan::stan_fit<model_sigfit_fitex_nmf_namespace::model_sigfit_fitex_nmf, boost::random::ecuyer1988> ::constrained_param_names)
    ;
}
