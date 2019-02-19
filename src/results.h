/**
 * Storage class for results.
 *
 * Author: Jiri Hon <jiri.hon@gmail.com>
 * Date: 2016/03/30
 * Package: pqsfinder
 */

#ifndef RESULTS_HEADER
#define RESULTS_HEADER

#include <Rcpp.h>
#include <cstdlib>
#include "features.h"
#include "opts.h"
#include "scores_buffer.h"

using namespace Rcpp;
using namespace std;


class results {
public:
  struct item_t {
    string::const_iterator start;
    int len;
    int score;
    int nt;
    int nb;
    int nm;
    int rl1;
    int rl2;
    int rl3;
    int ll1;
    int ll2;
    int ll3;
  };
  vector<results::item_t> items;
  
  int *density = NULL;
  int *max_scores = NULL;

  size_t seq_len = 0;
  string::const_iterator ref;
  
  scores_buffer scores;
  
  results(const size_t seq_len, const string::const_iterator ref, const opts_t &opts) :
    seq_len(seq_len), ref(ref), scores(opts.max_len, ref)
  {
    if (!opts.fast) {
      this->density = new int[seq_len];
      this->max_scores = new int[seq_len];
      
      for (size_t i = 0; i < seq_len; ++i) {
        this->density[i] = 0;
        this->max_scores[i] = 0;
      }
    }
  }
  ~results() {
    if (this->density != NULL)
      delete [] this->density;
    if (this->max_scores != NULL)
      delete [] this->max_scores;
  }
  inline void save_pqs(
      const int score, const string::const_iterator &s,
      const string::const_iterator &e, features_t &f)
  {
    results::item_t item;
    
    item.start = s;
    item.len = e - s;
    item.score = score;
    item.nt = f.nt;
    item.nb = f.nb;
    item.nm = f.nm;
    item.rl1 = f.rl1;
    item.rl2 = f.rl2;
    item.rl3 = f.rl3;
    item.ll1 = f.ll1;
    item.ll2 = f.ll2;
    item.ll3 = f.ll3;
    
    this->items.push_back(item);
  }
  inline void print() const {
    Rcout << "Results" << endl;
    for (unsigned i = 0; i < this->items.size(); i++) {
      Rcout << "PQS[" << i << "]: " << this->items[i].start - this->ref + 1 
            << "-" << this->items[i].start + this->items[i].len - this->ref << " "
            << string(this->items[i].start, this->items[i].start + this->items[i].len)
            << " " << this->items[i].score << endl;
    }
  }
  void sort_items() {
    sort(this->items.begin(), this->items.end(), results::compare_by_start);
  }
  void clear_max_scores() {
    if (this->max_scores != NULL) {
      for (size_t i = 0; i < this->seq_len; ++i) {
        this->max_scores[i] = 0;
      }
    }
  }
  static bool compare_by_start(const results::item_t &a, const results::item_t &b) {
    return a.start < b.start;
  }
};

#endif // RESULTS_HEADER
