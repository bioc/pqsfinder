/**
 * Algorithm options
 *
 * Author: Jiri Hon <jiri.hon@gmail.com>
 * Date: 2018/02/11
 * Package: pqsfinder
 */

#ifndef OPTS_HEADER
#define OPTS_HEADER


// algorithm options
struct opts_t {
  int max_len;
  int min_score;
  int run_min_len;
  int run_max_len;
  int loop_min_len;
  int loop_max_len;
  int check_int_period;
  size_t threads;
  size_t chunk_size;
  bool verbose;
  bool overlapping;
  bool use_re;
  bool debug;
  bool use_default_scoring;
  bool fast;
};


#endif // OPTS_HEADER
