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
  vector<int> loop_penalties;
  vector<int> bulge_penalties;
  
  scoring() {}
  ~scoring() {
    if (custom_scoring_fn != NULL) {
      delete custom_scoring_fn;
    }
  }
  void init_loop_penalties(int loop_max_len) {
    int len_sum_max = loop_max_len * 3;
    
    for (int i = 0; i < len_sum_max + 1; ++i) {
      double mean = ((double) i) / 3.0;
      loop_penalties.push_back((int) round(this->loop_mean_factor * pow(mean, this->loop_mean_exponent)));
    }
  }
  void init_bulge_penalties(int run_max_len) {
    // max_bulge_len = run_max_len - 2, two side guanines are obligatory
    
    for (int i = 0; i < run_max_len - 1; ++i) {
      bulge_penalties.push_back((int) round(this->bulge_len_factor * pow(i, this->bulge_len_exponent)));
    }
  }
  void set_custom_scoring_fn(SEXP custom_scoring_fn) {
    this->custom_scoring_fn = new Function(custom_scoring_fn);
  }
};


#endif // SCORING_HEADER
