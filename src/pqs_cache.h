/**
 * Cache for low complexity regions. It is usefull just for dealing with almost
 * G-complete sequence regions, which does not usually occur in human genome.
 *
 * Author: Jiri Hon <jiri.hon@gmail.com>
 * Date: 2016/03/30
 * Package: pqsfinder
 */

#ifndef PQS_CACHE_HEADER
#define PQS_CACHE_HEADER

#include <Rcpp.h>
#include <cstdlib>

using namespace Rcpp;
using namespace std;

class pqs_cache {
public:
  class entry {
  private:
    void init(const int max_len) {
      this->density = (int *)calloc(this->max_len, sizeof(int));
      if (this->density == NULL)
        stop("Unable to allocate memory for cache density vector.");
      this->score_dist = (int *)calloc(this->max_len, sizeof(int));
      if (this->score_dist == NULL)
        stop("Unable to allocate memory for cache score distribution vector.");
    }
  public:
    int *density;
    int *score_dist;
    int score;
    int len;
    const int max_len;
    entry(const int max_len) : score(0), len(0), max_len(max_len) {
      init(max_len);
    }
    entry(const entry &obj) : score(obj.score), len(obj.len), max_len(obj.max_len) {
      init(obj.max_len);
      memcpy(this->density, obj.density, this->max_len);
      memcpy(this->score_dist, obj.score_dist, this->max_len);
    }
    ~entry() {
      if (this->density != NULL)
        free(this->density);
      if (this->score_dist != NULL)
        free(this->score_dist);
    }
  };
  typedef map<string, pqs_cache::entry> cache_map;

  static const int use_treshold = 10000;
  cache_map table;
  int max_len;

  pqs_cache(const int max_len) : max_len(max_len) {}

  inline entry *get(const string::const_iterator &s, const string::const_iterator &e) {
    cache_map::iterator it = this->table.find(string(s, e));
    if (it == this->table.end())
      return NULL;
    else
      return &((*it).second);
  }
  inline void put(const string::const_iterator &s, const string::const_iterator &e,
                  const pqs_cache::entry &entry) {
    string seq = string(s, e);
    cache_map::iterator it = this->table.find(seq);
    if (it == this->table.end())
      this->table.insert(cache_map::value_type(seq, entry));
  }
};

#endif
