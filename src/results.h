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

using namespace Rcpp;
using namespace std;

class results {
public:
  struct item_t {
    int start;
    int len;
    int score;
    string strand;
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
  
  int *density;
  int *max_scores;

  const int min_score;
  const int seq_len;

  results(const int seq_len, const int min_score) :
    min_score(min_score), seq_len(seq_len)
  {
    this->density = (int *)calloc(seq_len, sizeof(int));
    if (this->density == NULL)
      throw runtime_error("Unable to allocate memory for results density vector.");
    this->max_scores = (int *)calloc(seq_len, sizeof(int));
    if (this->max_scores == NULL)
      throw runtime_error("Unable to allocate memory for results score distribution vector.");
  }
  ~results() {
    if (this->density != NULL)
      free(this->density);
    if (this->max_scores != NULL)
      free(this->max_scores);
  }
  inline void save_pqs(
      const int score, const string::const_iterator &s,
      const string::const_iterator &e, features_t &f,
      const string::const_iterator &ref, const string &strand)
  {
    if (score >= this->min_score) {
      results::item_t item;
      
      if (strand == "+") {
        item.start = s - ref + 1; // R indexing starts at 1
      } else {
        item.start = this->seq_len - (e - ref) + 1;
      }
      item.len = e - s;
      item.score = score;
      item.strand = strand;
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
  }
  inline void save_density_and_max_scores(
      const string::const_iterator &s, const string::const_iterator &ref,
      const string &strand, const int *density, const int *max_scores, const int max_len)
  {
    int offset, k_limit;
    if (strand == "+") {
      offset = s - ref;
      k_limit = min(max_len, this->seq_len - offset);
      for (int k = 0; k < k_limit; ++k) {
        int i = offset + k;
        this->density[i] += density[k];
        this->max_scores[i] = max(this->max_scores[i], max_scores[k]);
      }
    }
    else {
      offset = (this->seq_len - 1) - (s - ref);
      k_limit = min(max_len, offset + 1);
      for (int k = 0; k < k_limit; ++k) {
        int i = offset - k;
        this->density[i] += density[k];
        this->max_scores[i] = max(this->max_scores[i], max_scores[k]);
      }
    }
  }
  inline void print(const string::const_iterator &ref) const {
    Rcout << "Results" << endl;
    for (unsigned i = 0; i < this->items.size(); i++) {
      Rcout << "PQS[" << i << "]: " << this->items[i].start << " "
            << string(ref + this->items[i].start, ref + this->items[i].start + this->items[i].len)
            << " " << this->items[i].score << endl;
    }
  }
};

#endif // RESULTS_HEADER
