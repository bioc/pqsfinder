/**
 * Scoring options
 *
 * Author: Jiri Hon <jiri.hon@gmail.com>
 * Date: 2018/02/11
 * Package: pqsfinder
 */

#ifndef SCORING_HEADER
#define SCORING_HEADER

#include <Rcpp.h>

using namespace Rcpp;
using namespace std;


class scoring {
public:
  int tetrad_bonus;
  int bulge_penalty;
  double bulge_len_factor;
  double bulge_len_exponent;
  int mismatch_penalty;
  double loop_mean_factor;
  double loop_mean_exponent;
  int max_bulges;
  int max_mimatches;
  int max_defects;
  Function *custom_scoring_fn = NULL;
  int *loop_penalties = NULL;
  int *bulge_penalties = NULL;
  
  scoring() {}
  ~scoring() {
    if (custom_scoring_fn != NULL) {
      delete custom_scoring_fn;
    }
    if (loop_penalties != NULL) {
      delete [] loop_penalties;
    }
    if (bulge_penalties != NULL) {
      delete [] bulge_penalties;
    }
  }
  void init_loop_penalties(int loop_max_len) {
    int len_sum_max = loop_max_len * 3;
    loop_penalties = new int[len_sum_max + 1]; // to store len_sum_max + 1 items.
    
    for (int i = 0; i < len_sum_max + 1; ++i) {
      double mean = ((double) i) / 3.0;
      loop_penalties[i] = (int) round(this->loop_mean_factor * pow(mean, this->loop_mean_exponent));
    }
  }
  void init_bulge_penalties(int run_max_len) {
    bulge_penalties = new int[run_max_len - 1];
    // max_bulge_len = run_max_len - 2, two side guanines are obligatory
    
    for (int i = 0; i < run_max_len - 1; ++i) {
      bulge_penalties[i] = (int) round(this->bulge_len_factor * pow(i, this->bulge_len_exponent));
    }
  }
  void set_custom_scoring_fn(SEXP custom_scoring_fn) {
    this->custom_scoring_fn = new Function(custom_scoring_fn);
  }
};


#endif // SCORING_HEADER
