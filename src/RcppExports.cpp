// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <Rcpp.h>

using namespace Rcpp;

#ifdef RCPP_USE_GLOBAL_ROSTREAM
Rcpp::Rostream<true>&  Rcpp::Rcout = Rcpp::Rcpp_cout_get();
Rcpp::Rostream<false>& Rcpp::Rcerr = Rcpp::Rcpp_cerr_get();
#endif

// pqsfinder
SEXP pqsfinder(SEXP subject, std::string strand, bool overlapping, int max_len, int min_score, int run_min_len, int run_max_len, int loop_min_len, int loop_max_len, int max_bulges, int max_mismatches, int max_defects, int tetrad_bonus, int mismatch_penalty, int bulge_penalty, double bulge_len_factor, double bulge_len_exponent, double loop_mean_factor, double loop_mean_exponent, std::string run_re, SEXP custom_scoring_fn, bool use_default_scoring, bool deep, bool verbose);
RcppExport SEXP _pqsfinder_pqsfinder(SEXP subjectSEXP, SEXP strandSEXP, SEXP overlappingSEXP, SEXP max_lenSEXP, SEXP min_scoreSEXP, SEXP run_min_lenSEXP, SEXP run_max_lenSEXP, SEXP loop_min_lenSEXP, SEXP loop_max_lenSEXP, SEXP max_bulgesSEXP, SEXP max_mismatchesSEXP, SEXP max_defectsSEXP, SEXP tetrad_bonusSEXP, SEXP mismatch_penaltySEXP, SEXP bulge_penaltySEXP, SEXP bulge_len_factorSEXP, SEXP bulge_len_exponentSEXP, SEXP loop_mean_factorSEXP, SEXP loop_mean_exponentSEXP, SEXP run_reSEXP, SEXP custom_scoring_fnSEXP, SEXP use_default_scoringSEXP, SEXP deepSEXP, SEXP verboseSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< SEXP >::type subject(subjectSEXP);
    Rcpp::traits::input_parameter< std::string >::type strand(strandSEXP);
    Rcpp::traits::input_parameter< bool >::type overlapping(overlappingSEXP);
    Rcpp::traits::input_parameter< int >::type max_len(max_lenSEXP);
    Rcpp::traits::input_parameter< int >::type min_score(min_scoreSEXP);
    Rcpp::traits::input_parameter< int >::type run_min_len(run_min_lenSEXP);
    Rcpp::traits::input_parameter< int >::type run_max_len(run_max_lenSEXP);
    Rcpp::traits::input_parameter< int >::type loop_min_len(loop_min_lenSEXP);
    Rcpp::traits::input_parameter< int >::type loop_max_len(loop_max_lenSEXP);
    Rcpp::traits::input_parameter< int >::type max_bulges(max_bulgesSEXP);
    Rcpp::traits::input_parameter< int >::type max_mismatches(max_mismatchesSEXP);
    Rcpp::traits::input_parameter< int >::type max_defects(max_defectsSEXP);
    Rcpp::traits::input_parameter< int >::type tetrad_bonus(tetrad_bonusSEXP);
    Rcpp::traits::input_parameter< int >::type mismatch_penalty(mismatch_penaltySEXP);
    Rcpp::traits::input_parameter< int >::type bulge_penalty(bulge_penaltySEXP);
    Rcpp::traits::input_parameter< double >::type bulge_len_factor(bulge_len_factorSEXP);
    Rcpp::traits::input_parameter< double >::type bulge_len_exponent(bulge_len_exponentSEXP);
    Rcpp::traits::input_parameter< double >::type loop_mean_factor(loop_mean_factorSEXP);
    Rcpp::traits::input_parameter< double >::type loop_mean_exponent(loop_mean_exponentSEXP);
    Rcpp::traits::input_parameter< std::string >::type run_re(run_reSEXP);
    Rcpp::traits::input_parameter< SEXP >::type custom_scoring_fn(custom_scoring_fnSEXP);
    Rcpp::traits::input_parameter< bool >::type use_default_scoring(use_default_scoringSEXP);
    Rcpp::traits::input_parameter< bool >::type deep(deepSEXP);
    Rcpp::traits::input_parameter< bool >::type verbose(verboseSEXP);
    rcpp_result_gen = Rcpp::wrap(pqsfinder(subject, strand, overlapping, max_len, min_score, run_min_len, run_max_len, loop_min_len, loop_max_len, max_bulges, max_mismatches, max_defects, tetrad_bonus, mismatch_penalty, bulge_penalty, bulge_len_factor, bulge_len_exponent, loop_mean_factor, loop_mean_exponent, run_re, custom_scoring_fn, use_default_scoring, deep, verbose));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_pqsfinder_pqsfinder", (DL_FUNC) &_pqsfinder_pqsfinder, 24},
    {NULL, NULL, 0}
};

RcppExport void R_init_pqsfinder(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
