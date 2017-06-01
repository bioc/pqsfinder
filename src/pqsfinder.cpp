/**
 * Implementation of PQS search algorithm.
 *
 * Author: Jiri Hon <jiri.hon@gmail.com>
 * Date: 2016/01/17
 * Package: pqsfinder
 */

#include <Rcpp.h>
#include <string>
#include <climits>
#include <cstring>
#include <cstdio>
#include <cstdlib>
#include <algorithm>
#include <chrono>
#include <boost/regex.hpp>
#ifdef _GLIBCXX_DEBUG
#include <google/profiler.h>
#endif
#include "results.h"
#include "pqs_storage.h"
#include "pqs_cache.h"
#include "features.h"

using namespace Rcpp;
using namespace std;


/*
 * Extract from C++ reference documentation to string iterator semantics:
 *
 * string::end method returns an iterator pointing to the __past-the-end__ character
 * of the string! The past-the-end character is a theoretical character that would
 * follow the last character in the string. It __shall not be dereferenced.__
 *
 * The reason for that is because the ranges used by functions of the standard library
 * __do not include__ the element pointed by their closing iterator.
 */


// Implementation constants
static const int RUN_CNT = 4;


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
  Function *custom_scoring_fn;

  scoring() {
    custom_scoring_fn = NULL;
  }
  ~scoring() {
    if (custom_scoring_fn != NULL)
      delete custom_scoring_fn;
  }
  void set_custom_scoring_fn(SEXP custom_scoring_fn) {
    this->custom_scoring_fn = new Function(custom_scoring_fn);
  }
};

typedef struct flags {
  bool use_cache;
  bool use_re;
  bool use_prof;
  bool verbose;
  bool debug;
  bool use_default_scoring;
} flags_t;

typedef struct opts {
  bool overlapping;
  int max_len;
  int min_score;
  int run_min_len;
  int run_max_len;
  int loop_min_len;
  int loop_max_len;
  int check_int_period;
} opts_t;

// Class representing one run from quadruplex
class run_match {
public:
  string::const_iterator first;
  string::const_iterator second;
  int length() const {
    return second - first;
  };
};


/**
 * Print quadruplex summary
 *
 * @param m Quadruplex runs
 * @param score Score
 * @param ref Reference point, typically start of sequence
 */
inline void print_pqs(const run_match m[], int score, const string::const_iterator ref)
{
  Rcout << m[0].first - ref + 1 << "-" << m[3].second - m[0].first << " " << "[" << string(m[0].first, m[0].second) << "]";
  for (int i = 1; i < RUN_CNT; i++)
    Rcout << string(m[i-1].second, m[i].first) << "[" << string(m[i].first, m[i].second) << "]";
  Rcout << " " << score << endl;
}


/**
 * Count number of G's in G-run.
 * This way of counting leads to only one bulge allowed in G-run.
 *
 * @param m G-run match
 * @return Number of guanines in G-run
 */
inline int count_g_num(const run_match &m) {
  string::const_iterator s = m.first, e = m.second;
  int cnt = 0;
  while (*s == 'G' && s < e) { ++s; ++cnt; }
  --e; // <e> points to past-the-end character and as such should not be dereferenced
  while (*e == 'G' && e > s) { --e; ++cnt; }
  return cnt;
}


/**
 * R interface to test count_g_num function
 *
 * @param seq G-run sequence
 * @return Number of canonical Gs
 */
// ![[Rcpp::export]]
void count_g(std::string seq) {
  run_match m;
  m.first = seq.begin();
  m.second = seq.end();
  int cnt = count_g_num(m);
  Rcout << cnt << endl;
}


/**
 * Score run defects relative to the reference G-run.
 * 
 * @param pi Index of the reference G-run.
 * @param w Run widths.
 * @param g G contents.
 * @param l Run lenghts.
 * @param f PQS features.
 * @param sc Scoring options.
 * @return Scoring for G-runs.
 */
inline int score_run_defects(
    const int pi,
    const int w[],
    const int g[],
    features_t &f,
    const scoring &sc)
{
  int mismatches = 0, bulges = 0, perfects = 0;
  int score = 0;
  
  for (int i = 0; i < RUN_CNT; ++i) {
    if (w[i] == w[pi] && g[i] == g[pi]) {
      ++perfects;
    } else if ((w[i] == w[pi] && g[i] == g[pi] - 1)) {
      ++mismatches;
    } else if (w[i] > w[pi] && g[i] >= g[pi]) {
      ++bulges;
      score = score - (int) round(sc.bulge_len_factor * pow(w[i] - w[pi], sc.bulge_len_exponent));
    } else {
      return 0;
    }
  }
  if (mismatches <= sc.max_mimatches &&
      bulges <= sc.max_bulges &&
      mismatches + bulges <= sc.max_defects)
  {
    score = score + (w[pi] - 1) * sc.tetrad_bonus
            - mismatches * sc.mismatch_penalty
            - bulges * sc.bulge_penalty;
    f.nt = w[pi];
    f.nb = bulges;
    f.nm = mismatches;
    return score;
  } else {
    return 0;
  }
}

/**
 * Score PQS.
 *
 * @param m Quadruples runs.
 * @param f PQS features.
 * @param sc Scoring table.
 * @param opts Algorithm options.
 * @return Quadruplex score.
 */
inline int score_pqs(
    run_match m[],
    features_t &f,
    const scoring &sc,
    const opts_t &opts)
{
  int w[RUN_CNT], g[RUN_CNT], l[RUN_CNT - 1];
  int min_pw, pi, score;
  double mean;
  
  l[0] = m[1].first - m[0].second;
  l[1] = m[2].first - m[1].second;
  l[2] = m[3].first - m[2].second;
  
  // check if no more than one loop has zero length
  if (opts.loop_min_len == 0 &&
      ( (l[0] == 0 && l[1] == 0) ||
        (l[0] == 0 && l[2] == 0) ||
        (l[1] == 0 && l[2] == 0) ) ) {
    return 0;
  }

  w[0] = m[0].length();
  w[1] = m[1].length();
  w[2] = m[2].length();
  w[3] = m[3].length();
  
  g[0] = count_g_num(m[0]);
  g[1] = count_g_num(m[1]);
  g[2] = count_g_num(m[2]);
  g[3] = count_g_num(m[3]);
  
  min_pw = INT_MAX;
  pi = -1;
  for (int i = 0; i < RUN_CNT; ++i) {
    if (w[i] < min_pw && g[i] == w[i]) {
      min_pw = w[i];
      pi = i;
    }
  }
  if (pi < 0) {
    return 0;
  }

  score = score_run_defects(pi, w, g, f, sc);
  if (score <= 0) {
    return 0;
  }
  // update reported loop lengths
  f.ll1 = l[0];
  f.ll2 = l[1];
  f.ll3 = l[2];
  
  mean = (double) (l[0] + l[1] + l[2]) / 3.0;
  
  return max(score - (int) round(sc.loop_mean_factor * pow(mean, sc.loop_mean_exponent)), 0);
}


/**
 * Check user scoring function
 *
 * @param score Quadruplex score
 * @param m Quadruples runs
 * @param sc Scoring table
 * @example user function in R
   my_fn <- function(
      subject, score, start, width, loop_1, run_2, loop_2,
      run_3, loop_3, run_4)
   {
      len <- loop_1 - start
      if (len == loop_2 - run_2 && len == loop_3 - run_3 &&
          len == start + width - run_4)
      return(200)
   }
 */
inline void check_custom_scoring_fn(
    int &score, const run_match m[], const scoring &sc, SEXP subject,
    const string::const_iterator ref)
{
  int start, width, loop_1, run_2, loop_2, run_3, loop_3, run_4;

  start = m[0].first - ref + 1;
  width = m[3].second - m[0].first;
  loop_1 = m[0].second - ref + 1;
  run_2 = m[1].first - ref + 1;
  loop_2 = m[1].second - ref + 1;
  run_3 = m[2].first - ref + 1;
  loop_3 = m[2].second - ref + 1;
  run_4 = m[3].first - ref + 1;

  score = as<int>((*sc.custom_scoring_fn)(
    subject, score, start, width,
    loop_1, run_2, loop_2, run_3, loop_3, run_4));
}


/**
 * Run Boost regex engine.
 * 
 * @param start Start position.
 * @param end End position.
 * @param boost_m Boost output array.
 * @param run_re_c Regular expression.
 * @return Status.
 */
inline bool run_regex_search(
    const string::const_iterator &start,
    const string::const_iterator &end,
    boost::smatch &boost_m,
    const boost::regex &run_re_c)
{
  try {
    return boost::regex_search(start, end, boost_m, run_re_c, boost::match_default);
  } catch (bad_alloc &ba) {
    throw runtime_error(string("Regexp engine failed with exception: ") + ba.what());
    return false;
  }
}

/**
 * Perform run search on particular sequence region
 *
 * @param s Start of region
 * @param e End of region (iterator pointing past-the-end)
 * @param m Match info structure
 * @param run_re_c Run regular expression
 * @param opts Algorithm options and limits
 * @param flags Algorithm flags
 * @return True on success, false otherwise
 */
inline bool find_run(
    const string::const_iterator &start,
    const string::const_iterator &end,
    run_match &m,
    const boost::regex &run_re_c,
    const opts_t &opts,
    const flags_t &flags)
{
  if (flags.use_re) {
    static boost::smatch boost_m;
    string::const_iterator s, e;
    bool status = run_regex_search(start, end, boost_m, run_re_c);
    
    if (status) {
      if (boost_m[0].second - boost_m[0].first > opts.run_max_len) {
        // too long run found, search again with limited end
        s = boost_m[0].first;
        e = min(s + opts.run_max_len, end);
        status = run_regex_search(s, e, boost_m, run_re_c);
        if (!status) {
          return false;
        }
      }
      if (boost_m[0].second - boost_m[0].first < opts.run_min_len) {
        return false;
      }
      m.first = boost_m[0].first;
      m.second = boost_m[0].second;
      return true;
    } else {
      return false;
    }
  } else {
    string::const_iterator s = start, e;
    
    while (*s != 'G' && s < end) ++s;
    e = min(s + opts.run_max_len, end);
    --e; // <e> points to past-the-end character and as such should not be dereferenced
    while (*e != 'G' && e > s) --e;

    if (e - s + 1 < opts.run_min_len) {
      // definitely too short to be a proper run
      return false;
    }
    m.first = s;
    m.second = ++e; // correction to point on the past-the-end character
    return true;
  }
}


/**
 * Debug start and end coordinates during the iterative process.
 * 
 * @param name Label
 * @param i Start index
 * @param s Start of region
 * @param e End of region
 * @param ref Reference point, typically start of sequence
 */
void debug_s_e(
    const char *name,
    int i,
    string::const_iterator &s,
    string::const_iterator &e,
    const string::const_iterator &ref) {
  
  int s_i = s - ref + 1;
  int e_i = e - ref;

  if (s_i == 6 || s_i == 7) {
    Rprintf("[%d] %s: %d %d\n", i, name, s_i, e_i);
  }
}


/**
 * Recursively idetify 4 consecutive runs making quadruplex
 *
 * @param subject DNAString object
 * @param strand Strand specification
 * @param i Odinal number of quadruplex run
 * @param start Start position for the current run
 * @param end Limit end position for the current run
 * @param m Array of run matches
 * @param run_re_c Compiled run regular expression
 * @param opts Algorithm options
 * @param sc Scoring options
 * @param ref Reference point, typically start of sequence
 * @param len Total sequence length
 * @param pqs_start Start of the first G-run
 * @param pqs_storage Storage object
 * @param ctable Cache table
 * @param cache_entry Candidate cache entry
 * @param pqs_cnt PQS counter
 * @param res Object for output PQS
 * @param zero_loop Flag, if PQS has zero-length loop
 */
void find_all_runs(
    SEXP subject,
    const string &strand,
    const int i,
    string::const_iterator start,
    string::const_iterator end,
    run_match m[],
    const boost::regex &run_re_c,
    const opts_t &opts,
    const flags_t &flags,
    const scoring &sc,
    const string::const_iterator &ref,
    const size_t len,
    string::const_iterator &pqs_start,
    pqs_storage &pqs_storage,
    pqs_cache &ctable,
    pqs_cache::entry &cache_entry,
    int &pqs_cnt,
    results &res,
    bool zero_loop,
    chrono::system_clock::time_point s_time)
{
  string::const_iterator s, e, min_e;
  int score, loop_len;
  pqs_cache::entry *cache_hit;
  
  if (i > 0) {
    loop_len = start - m[i-1].second; 
    if (loop_len < opts.loop_min_len) {
      start = min(m[i-1].second + opts.loop_min_len, end); // skip too short loop
    } else if (flags.use_default_scoring && loop_len == 0 && zero_loop) {
      start = min(m[i-1].second + 1, end); // only one zero-length loop is allowed
    }
  }
  
  for (s = start; s < end; ++s)
  {
    if (i == 0)
    {// specific code for the first run matching
      if (flags.use_cache && cache_entry.density[0] > pqs_cache::use_treshold)
      {
        cache_hit = ctable.get(s, min(s + opts.max_len, end));

        if (cache_hit != NULL) {
          if (flags.debug)
            Rcout << "Cache hit: " << s - ref  << " " << string(s, s+cache_hit->len)
                  << " " << cache_hit->score << endl;

          res.save_density_and_max_scores(
            s, ref, strand, cache_hit->density, cache_hit->max_scores, opts.max_len);

          pqs_storage.insert_pqs(cache_hit->score, s, s + cache_hit->len, cache_hit->f, res, ref, strand);
          continue;
        }
      }
      // reset score of best PQS starting at current position
      cache_entry.score = 0;
      // reset density and score distribution
      for (int k = 0; k < opts.max_len; ++k) {
        cache_entry.density[k] = 0;
        cache_entry.max_scores[k] = 0;
      }
    }
    min_e = s + opts.run_min_len;

    for (e = end; e >= min_e && find_run(s, e, m[i], run_re_c, opts, flags); e--)
    {
      // update search bounds
      s = string::const_iterator(m[i].first);
      e = string::const_iterator(m[i].second);
      
      loop_len = s - m[i-1].second;
      if (loop_len == 0) {
        zero_loop = true;
      }
      if (i > 0 && loop_len > opts.loop_max_len)
        return; // skip too long loops

      if (i == 0)
        // enforce G4 total length limit to be relative to the first G-run start
        find_all_runs(
          subject, strand, i+1, e, min(s + opts.max_len, end), m, run_re_c,
          opts, flags, sc, ref, len, s, pqs_storage, ctable, cache_entry,
          pqs_cnt, res, zero_loop, s_time
        );
      else if (i < 3)
        find_all_runs(
          subject, strand, i+1, e, end, m, run_re_c, opts, flags, sc, ref, len,
          pqs_start, pqs_storage, ctable, cache_entry, pqs_cnt, res, zero_loop, s_time
        );
      else {
        /* Check user interrupt after reasonable amount of PQS identified to react
         * on important user signals. I.e. he might want to abort the computation. */
        if (++pqs_cnt == opts.check_int_period)
        {
          pqs_cnt = 0;
          checkUserInterrupt();
          if (!flags.verbose) {
            double percents = ceilf((m[0].first - ref)/(double)len*100);
            double e_seconds = chrono::duration_cast<std::chrono::seconds>(
              chrono::system_clock::now() - s_time).count();
            double r_seconds = (e_seconds / percents) * (100 - percents);
            int r_hours = r_seconds / 3600;
            r_seconds -= r_hours * 3600;
            int r_minutes = r_seconds / 60;
            r_seconds -= r_minutes * 60;
            
            char buffer[10];
            sprintf(buffer, "%02d:%02d:%02d", r_hours, r_minutes, (int) r_seconds);
            
            Rcout << "Search status: " << percents << "% ETTC " << string(buffer) << "\r" << flush;
          }
        }
        score = 0;
        features_t pqs_features;
        if (flags.use_default_scoring) {
          score = score_pqs(m, pqs_features, sc, opts);
        }
        if ((score || !flags.use_default_scoring) && sc.custom_scoring_fn != NULL) {
          check_custom_scoring_fn(score, m, sc, subject, ref);
        }
        if (score) {
          int pqs_len = m[3].second - m[0].first;
          
          for (int k = 0; k < pqs_len; ++k) {
            cache_entry.max_scores[k] = max(cache_entry.max_scores[k], score);
          }
          if (score >= opts.min_score) {
            // current PQS satisfied all constraints
            pqs_storage.insert_pqs(score, m[0].first, m[3].second, pqs_features, res, ref, strand);
            
            for (int k = 0; k < pqs_len; ++k)
              ++cache_entry.density[k];
            
            if (score > cache_entry.score ||
                (score == cache_entry.score && pqs_len < cache_entry.len)) {
              // update properties of caching candidate
              cache_entry.score = score;
              cache_entry.len = pqs_len;
              cache_entry.f = pqs_features;
            }
            if (flags.verbose)
              print_pqs(m, score, ref);
          }
        }
      }
    }
    if (i == 0) {
      if (flags.use_cache && cache_entry.density[0] > pqs_cache::use_treshold)
        ctable.put(s, min(s + opts.max_len, end), cache_entry);

      // add locally accumulated max scores to global max scores array
      res.save_density_and_max_scores(
        s, ref, strand, cache_entry.density, cache_entry.max_scores, opts.max_len);
    }
  }
}

/**
 * Select between overlapping and non-overlapping storage.
 * 
 * @param overlapping If report overlapping PQS.
 * @param ov Overlapping storage.
 * @param nov Non-overlapping storage.
 * @return Reference to storage interface.
 */
pqs_storage &select_pqs_storage(
    bool overlapping,
    pqs_storage_overlapping &ov,
    pqs_storage_non_overlapping &nov)
{
  if (overlapping)
    return ov;
  else
    return nov;
}

/**
 * Perform quadruplex search on given DNA sequence.
 *
 * @param subject DNAString object
 * @param seq DNA sequence
 * @param strand Strand specification
 * @param run_re_c Run regular expression
 * @param ctable PQS cache table
 * @param sc Scoring options
 * @param opts Algorihtm options
 * @param flags Algorithm flags
 * @param res Results object
 */
void pqs_search(
    SEXP subject,
    const string &seq,
    const string strand,
    const boost::regex &run_re_c,
    pqs_cache &ctable,
    const scoring &sc,
    const opts_t &opts,
    const flags_t &flags,
    results &res)
{
  run_match m[RUN_CNT];
  pqs_cache::entry cache_entry(opts.max_len);
  string::const_iterator pqs_start;
  int pqs_cnt = 0;
  pqs_storage_overlapping pqs_storage_ov(seq.begin());
  pqs_storage_non_overlapping pqs_storage_nov;
  pqs_storage &pqs_storage = select_pqs_storage(opts.overlapping, pqs_storage_ov, pqs_storage_nov);
  
  // Global sequence length is the only limit for the first G-run
  find_all_runs(
    subject, strand, 0, seq.begin(), seq.end(), m, run_re_c, opts, flags, sc,
    seq.begin(), seq.length(), pqs_start, pqs_storage, ctable,
    cache_entry, pqs_cnt, res, 0, chrono::system_clock::now()
  );
  pqs_storage.export_pqs(res, seq.begin(), strand);
}


//' Identify potential quadruplex forming sequences.
//'
//' Function for identification of all potential intramolecular quadruplex
//' patterns (PQS) in DNA sequence.
//' 
//' Use \code{\link{elementMetadata}} function to get extra PQS features
//' like number of tetrads (nt), bulges (nb), mismatches (nm) or loop lengths
//' (ll1, ll2, ll3).
//'
//' @param subject DNAString object.
//' @param strand Strand specification. Allowed values are "+", "-" or "*",
//'   where the last one represents both strands. Implicitly, the input
//'   DNAString object is assumed to encode the "+" strand.
//' @param overlapping If true, than all overlapping PQS will be reported.
//' @param max_len Maximal lenth of PQS.
//' @param min_score Minimal PQS score.
//' @param run_min_len Minimal length of quadruplex run.
//' @param run_max_len Maximal length of quadruplex run.
//' @param loop_min_len Minimal length of quadruplex loop. Unless the default scoring
//' system is disabled, at most one loop can have zero length.
//' @param loop_max_len Maxmimal length of quadruplex loop.
//' @param max_bulges Maximal number of runs with bulge.
//' @param max_mismatches Maximal number of runs with mismatch.
//' @param max_defects Maximum number of defects in total (\code{max_bulges +
//'   max_mismatches}).
//' @param tetrad_bonus Score bonus for one complete G tetrade.
//' @param mismatch_penalty Penalization for a mismatch in tetrad.
//' @param edge_mismatch_penalty Penalization for an mismatch in edge tetrad.
//' @param bulge_penalty Penalization for a bulge in quadruplex run.
//' @param bulge_len_factor Penalization factor for a bulge length.
//' @param bulge_len_exponent Exponent of bulge length.
//' @param loop_mean_factor Penalization factor of loop length mean.
//' @param loop_mean_exponent Exponent of loop length mean.
//' @param run_re Regular expression specifying one run of quadruplex.
//' @param custom_scoring_fn Custom quadruplex scoring function. It takes the
//'   following 10 arguments: \code{subject} - Input DNAString object,
//'   \code{score} - implicit PQS score, \code{start} - PQS start position,
//'   \code{width} - PQS width, \code{loop_1} - start pos. of loop #1,
//'   \code{run_2} - start pos. of run #2, \code{loop_2} - start pos. of loop
//'   #2, \code{run_3} - start pos. of run #3, \code{loop_3} - start pos. of
//'   loop #3, \code{run_4} - start pos. of run #4. Return value of the function
//'   has to be new score represented as a single integer value. Please note
//'   that if \code{use_default_scoring} is enabled, the custom scoring function
//'   is evaluated AFTER the default scoring system but ONLY IF the default
//'   scoring system resulted in non-zero score (for performance reasons). On
//'   the other hand, when \code{use_default_scoring} is disabled, custom
//'   scoring function is evaluated on every PQS.
//' @param use_default_scoring Enables default internal scoring system. This
//'   option is particularly useful in case you intend to radically change the
//'   default behavior and specify your own scoring function. By disabling the
//'   default scoring you will get a full control above the underlying detection
//'   algorithm.
//' @param verbose Enables detailed output. Turn it on if you want to see all
//'   possible PQS found at each positions and not just the best one. It is
//'   highly recommended to use this option for debugging custom quadruplex
//'   scoring function. Each PQS is reported on separate row in the following
//'   format: \code{start cnt pqs_sequence score}, where \code{start} is the PQS
//'   starting position, \code{pqs_sequence} shows the PQS sequence structure
//'   with each run surrounded by square brackets and \code{score} is the score
//'   assigned to the particular PQS by all applied scoring functions.
//' @return \code{\link{PQSViews}} object
//'
//' @examples
//' pv <- pqsfinder(DNAString("CCCCCCGGGTGGGTGGGTGGTAAAA"))
//' pv
//' elementMetadata(pv)
//'
// [[Rcpp::export]]
SEXP pqsfinder(
    SEXP subject,
    std::string strand = "*",
    bool overlapping = false,
    int max_len = 50,
    int min_score = 26,
    int run_min_len = 2,
    int run_max_len = 11,
    int loop_min_len = 0,
    int loop_max_len = 30,
    int max_bulges = 3,
    int max_mismatches = 3,
    int max_defects = 3,
    int tetrad_bonus = 40,
    int mismatch_penalty = 28,
    int bulge_penalty = 20,
    double bulge_len_factor = 0.2,
    double bulge_len_exponent = 1,
    double loop_mean_factor = 6.6,
    double loop_mean_exponent = 0.8,
    std::string run_re = "G{1,10}.{0,9}G{1,10}",
    SEXP custom_scoring_fn = R_NilValue,
    bool use_default_scoring = true,
    bool verbose = false)
{
  if (max_len < 1)
    throw invalid_argument("Maximal length of PQS has to be a positive value.");
  if (min_score < 1) {
    throw invalid_argument("Minimal PQS score has to be a positive value.");
  }

  if (run_min_len < 2)
    throw invalid_argument("Minimal PQS run length has to be greater than 2 or equal.");
  if (run_max_len < 2)
    throw invalid_argument("Maximal PQS run length has to be greater than 2 or equal.");
  if (run_min_len > run_max_len)
    throw invalid_argument("Minimal PQS run length can't be greater than the maximal PQS run length.");

  if (loop_min_len < 0)
    throw invalid_argument("Minimal PQS loop length has to be a non-negative value.");
  if (loop_max_len < 0)
    throw invalid_argument("Maximal PQS loop length has to be a non-negative value.");
  if (loop_min_len > loop_max_len)
    throw invalid_argument("Minimal PQS loop length can't be greater than the maximal PQS loop length.");

  if (max_bulges < 0 || max_bulges > 3)
    throw invalid_argument("Maximum number of runs with bulges has to be from the range 0-3.");
  if (max_mismatches < 0 || max_mismatches > 3)
    throw invalid_argument("Maximum number of runs with mismatches has to be from the range 0-3.");
  if (max_defects < 0 || max_defects > 3)
    throw invalid_argument("Maximum number of runs with defects (bulge or mismatch) has to be from the range 0-3.");
  if (strand != "+" && strand != "-" && strand != "*")
    throw invalid_argument("Strand specification must be +, - or *.");

  Function as_character("as.character");
  Function get_class("class");

  CharacterVector subject_class = as_character(get_class(subject));

  if (subject_class[0] != "DNAString")
    throw invalid_argument("Subject must be DNAString object.");

  flags_t flags;
  flags.use_cache = overlapping ? false : true;
  flags.use_re = false;
  flags.use_prof = false;
  flags.debug = false;
  flags.verbose = verbose;
  flags.use_default_scoring = use_default_scoring;

  if (run_re != "G{1,10}.{0,9}G{1,10}") {
    // User specified its own regexp, force to use regexp engine
    flags.use_re = true;
  }

  opts_t opts;
  opts.overlapping = overlapping;
  opts.max_len = max_len;
  opts.min_score = min_score;
  opts.loop_max_len = loop_max_len;
  opts.loop_min_len = loop_min_len;
  opts.run_max_len = run_max_len;
  opts.run_min_len = run_min_len;
  
  if (flags.use_re) {
    opts.check_int_period = 1e6;
  } else {
    opts.check_int_period = 1e7;
  }
  scoring sc;
  sc.tetrad_bonus = tetrad_bonus;
  sc.bulge_penalty = bulge_penalty;
  sc.bulge_len_factor = bulge_len_factor;
  sc.bulge_len_exponent = bulge_len_exponent;
  sc.mismatch_penalty = mismatch_penalty;
  sc.loop_mean_factor = loop_mean_factor;
  sc.loop_mean_exponent = loop_mean_exponent;
  sc.max_bulges = max_bulges;
  sc.max_mimatches = max_mismatches;
  sc.max_defects = max_defects;

  string seq = as<string>(as_character(subject));
  Function reverseComplement("reverseComplement");
  SEXP subject_rc = reverseComplement(subject);
  string seq_rc = as<string>(as_character(subject_rc));

  results res(seq.length(), opts.min_score);
  pqs_cache ctable(opts.max_len);
  boost::regex run_re_c(run_re);

  if (flags.debug) {
    Rcout << "G-run regexp: " << run_re << endl;
    Rcout << "Use cache: " << flags.use_cache << endl;
    Rcout << "Use regexp engine: " << flags.use_re << endl;
    Rcout << "Input sequence length: " << seq.length() << endl;
    Rcout << "Use user fn: " << (custom_scoring_fn != R_NilValue) << endl;
  }

  if (custom_scoring_fn != R_NilValue) {
    sc.set_custom_scoring_fn(custom_scoring_fn);
    if (flags.debug) {
      Rcout << "User function: " << endl;
      Function show("show");
      show(custom_scoring_fn);
    }
  }

  #ifdef _GLIBCXX_DEBUG
  if (flags.use_prof)
    ProfilerStart("profiling.log");
  #endif

  if (strand == "+" || strand == "*") {
    Rcout << "Searching on sense strand..." << endl;
    pqs_search(subject, seq, "+", run_re_c, ctable, sc, opts, flags, res);
  }
  if (strand == "-" || strand == "*") {
    Rcout << "Searching on antisense strand..." << endl;
    pqs_search(subject_rc, seq_rc, "-", run_re_c, ctable, sc, opts, flags, res);
  }

  #ifdef _GLIBCXX_DEBUG
  if (flags.use_prof)
    ProfilerStop();
  #endif

  IntegerVector res_start(res.start.begin(), res.start.end());
  IntegerVector res_width(res.len.begin(), res.len.end());
  IntegerVector res_score(res.score.begin(), res.score.end());
  CharacterVector res_strand(res.strand.begin(), res.strand.end());
  IntegerVector res_nt(res.nt.begin(), res.nt.end());
  IntegerVector res_nb(res.nb.begin(), res.nb.end());
  IntegerVector res_nm(res.nm.begin(), res.nm.end());
  IntegerVector res_ll1(res.ll1.begin(), res.ll1.end());
  IntegerVector res_ll2(res.ll2.begin(), res.ll2.end());
  IntegerVector res_ll3(res.ll3.begin(), res.ll3.end());

  IntegerVector res_density(seq.length());
  IntegerVector res_max_scores(seq.length());
  for (unsigned i = 0; i < seq.length(); ++i) {
    res_density[i] = res.density[i];
    res_max_scores[i] = res.max_scores[i];
  }
  Function pqsviews("PQSViews");
  return pqsviews(
    subject, res_start, res_width, res_strand, res_score,
    res_density, res_max_scores,
    res_nt, res_nb, res_nm, res_ll1, res_ll2, res_ll3);
}
