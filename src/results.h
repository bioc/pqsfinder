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
  vector<int> start;
  vector<int> len;
  vector<int> score;
  vector<string> strand;
  vector<int> nt;
  vector<int> nb;
  vector<int> nm;
  vector<int> rl1;
  vector<int> rl2;
  vector<int> rl3;
  vector<int> ll1;
  vector<int> ll2;
  vector<int> ll3;
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
      if (strand == "+")
        this->start.push_back(s - ref + 1); // R indexing starts at 1
      else
        this->start.push_back(this->seq_len - (e - ref) + 1);

      this->len.push_back(e - s);
      this->score.push_back(score);
      this->strand.push_back(strand);
      this->nt.push_back(f.nt);
      this->nb.push_back(f.nb);
      this->nm.push_back(f.nm);
      this->rl1.push_back(f.rl1);
      this->rl2.push_back(f.rl2);
      this->rl3.push_back(f.rl3);
      this->ll1.push_back(f.ll1);
      this->ll2.push_back(f.ll2);
      this->ll3.push_back(f.ll3);
    }
  }
  inline void save_pqs(const results &res, const int i)
  {
    this->start.push_back(res.start[i]);
    this->len.push_back(res.len[i]);
    this->score.push_back(res.score[i]);
    this->strand.push_back(res.strand[i]);
    this->nt.push_back(res.nt[i]);
    this->nb.push_back(res.nb[i]);
    this->nm.push_back(res.nm[i]);
    this->rl1.push_back(res.rl1[i]);
    this->rl2.push_back(res.rl2[i]);
    this->rl3.push_back(res.rl3[i]);
    this->ll1.push_back(res.ll1[i]);
    this->ll2.push_back(res.ll2[i]);
    this->ll3.push_back(res.ll3[i]);
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
    for (unsigned i = 0; i < this->start.size(); i++) {
      Rcout << "PQS[" << i << "]: " << this->start[i] << " "
            << string(ref + this->start[i], ref + this->start[i] + this->len[i])
            << " " << this->score[i] << endl;
    }
  }
};

#endif // RESULTS_HEADER
