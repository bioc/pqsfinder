/**
 * pqsfinder: an exhaustive and imperfection-tolerant search tool for potential
 * quadruplex-forming sequences (PQS)
 *
 * Author: Jiri Hon <jiri.hon@gmail.com>,
 *         Dominika Labudova <d.labudova@gmail.com>
 * Date: 2016/01/17
 * Package: pqsfinder
 */

#define GPERF_ENABLED

#include <Rcpp.h>
#include <string>
#include <climits>
#include <cstring>
#include <cstdio>
#include <cstdlib>
#include <algorithm>
#include <chrono>
#include <thread>
#include <boost/regex.hpp>
#include <boost/thread.hpp>
#ifdef GPERF_ENABLED
#include <gperftools/profiler.h>
#endif
#include "overscored.h"
#include "pqsfinder.h"

using namespace Rcpp;
using namespace std;


/*
 * Excerpt from C++ reference documentation on string iterator semantics:
 *
 * string::end method returns an iterator pointing to the __past-the-end__ character
 * of the string! The past-the-end character is a theoretical character that would
 * follow the last character in the string. It __shall not be dereferenced.__
 *
 * The reason for that is because the ranges used by functions of the standard library
 * __do not include__ the element pointed by their closing iterator.
 */

/*
 * TODO:
 */



// search progress description
struct search_progress_t {
  int seconds;
  int minutes;
  int hours;
  double percents;
};


// sequence chunk
struct seq_chunk_t {
  string::const_iterator s;
  string::const_iterator e;
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
  Rcout << m[0].first - ref + 1 << "-" << m[3].second - ref << " " << "[" << string(m[0].first, m[0].second) << "]";
  for (int i = 1; i < RUN_CNT; i++)
    Rcout << string(m[i-1].second, m[i].first) << "[" << string(m[i].first, m[i].second) << "]";
  Rcout << " " << score << endl;
}


/**
 * Print partial quadruplex
 * 
 */
inline void print_partial_pqs(const run_match m[], int i, const string::const_iterator ref) {
  Rcout << m[0].first - ref + 1 << "-" << m[i].second - ref << " " << "[" << string(m[0].first, m[0].second) << "]";
  for (int k = 1; k <= i; k++)
    Rcout << string(m[k-1].second, m[k].first) << "[" << string(m[k].first, m[k].second) << "]";
  Rcout << endl;
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
      // score = score - (int) round(sc.bulge_len_factor * pow(w[i] - w[pi], sc.bulge_len_exponent));
      score = score - sc.bulge_penalties[w[i] - w[pi]];
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
  // double mean;
  
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
  
  f.rl1 = w[0];
  f.rl2 = w[1];
  f.rl3 = w[2];
  
  // update reported loop lengths
  f.ll1 = l[0];
  f.ll2 = l[1];
  f.ll3 = l[2];
  
  // mean = (double) (l[0] + l[1] + l[2]) / 3.0;
  
  // return max(score - (int) round(sc.loop_mean_factor * pow(mean, sc.loop_mean_exponent)), 0);
  return max(score - sc.loop_penalties[l[0] + l[1] + l[2]], 0);
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
 * @return True on success, false otherwise
 */
inline bool find_run(
    const string::const_iterator &start,
    const string::const_iterator &end,
    run_match &m,
    const boost::regex &run_re_c,
    const opts_t &opts)
{
  string::const_iterator s = start, e;
  
  if (opts.use_re) {
    static boost::smatch boost_m;
    bool status = false;
    
    while (s < end) {
      status = run_regex_search(s, end, boost_m, run_re_c);
      if (!status) {
        break;
      }
      if (boost_m[0].second - boost_m[0].first > opts.run_max_len) {
        s = boost_m[0].first;
        e = min(s + opts.run_max_len, end);
        status = run_regex_search(s, e, boost_m, run_re_c);
        if (status) {
          break;
        } else {
          ++s;
        }
      } else {
        break;
      }
    }
    if (status) {
      if (boost_m[0].second - boost_m[0].first < opts.run_min_len) {
        return false;
      }
      m.first = boost_m[0].first;
      m.second = boost_m[0].second;
      
      if (opts.fast) {
        m.g_count = count_g_num(m);
      }
      return true;
    } else {
      return false;
    }
  } else {
    while (s < end) {
      while (*s != 'G' && s < end) ++s;
      e = min(s + opts.run_max_len, end);
      --e; // <e> points to past-the-end character and as such should not be dereferenced
      while (*e != 'G' && e > s) --e;
      
      if (e - s + 1 >= opts.run_min_len) {
        break;
      } else {
        ++s;
      }
    }
    if (e - s + 1 < opts.run_min_len) {
      // definitely too short to be a proper run
      return false;
    }
    m.first = s;
    m.second = ++e; // correction to point on the past-the-end character
    
    if (opts.fast) {
      m.g_count = count_g_num(m);
    }
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
 * Compute search progress
 * 
 * @param s Current position in the sequence
 * @param ref The beginning of the sequence
 * @param len Length of the sequence
 * @param s_time Starting time
 */
inline search_progress_t search_progress(
    const string::const_iterator s,
    const string::const_iterator ref,
    const size_t len,
    const chrono::system_clock::time_point s_time)
{
  double percents = ceilf((s - ref)/(double)len*100);
  double e_seconds = chrono::duration_cast<std::chrono::seconds>(
    chrono::system_clock::now() - s_time).count();
  double r_seconds = (e_seconds / percents) * (100 - percents);
  int r_hours = r_seconds / 3600;
  r_seconds -= r_hours * 3600;
  int r_minutes = r_seconds / 60;
  r_seconds -= r_minutes * 60;
  
  search_progress_t sp;
  sp.seconds = r_seconds;
  sp.minutes = r_minutes;
  sp.hours = r_hours;
  sp.percents = percents;
  return sp;
}


/**
 * Recursively idetify 4 consecutive runs making quadruplex
 *
 * @param subject DNAString object
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
 * @param vec_cache Candidate cache entry
 * @param int_cnt PQS counter
 * @param res Object for output PQS
 * @param zero_loop Flag, if PQS has zero-length loop
 * @param s_time Starting time
 */
void find_all_runs(
    const SEXP subject,
    const int i,
    string::const_iterator start,
    string::const_iterator end,
    run_match m[],
    const boost::regex &run_re_c,
    const opts_t &opts,
    const scoring &sc,
    const string::const_iterator &ref,
    const size_t len,
    storage &pqs_storage,
    int &int_cnt,
    results &res,
    bool zero_loop,
    const chrono::system_clock::time_point s_time,
    int tetrad_count,
    int defect_count,
    int &fn_call_count)
{
  string::const_iterator s, e, min_e;
  int score, loop_len;
  bool found_any;
  int next_tetrad_count = INT_MAX;
  int next_defect_count = 0;
  int max_score = 0;
  
  ++fn_call_count;
  ++int_cnt;
  
  /* Check user interrupt after reasonable amount of runs matched
   * on important user signals. E.g. the user might want to abort the computation. */
  if (opts.threads == 1 && int_cnt > opts.check_int_period)
  {
    int_cnt = 0;
    checkUserInterrupt();
    
    if (!opts.verbose) {
      search_progress_t sp = search_progress(start, ref, len, s_time);
      char buffer[10];
      sprintf(buffer, "%02d:%02d:%02d", sp.hours, sp.minutes, sp.seconds);
      Rcout << "Search status: " << sp.percents << "% ETTC " << string(buffer) << "\r" << flush;
    }
  }

  if (i > 0) {
    loop_len = start - m[i-1].second;
    if (loop_len < opts.loop_min_len) {
      start = min(m[i-1].second + opts.loop_min_len, end); // skip too short loop
    } else if (opts.use_default_scoring && loop_len == 0 && zero_loop) {
      start = min(m[i-1].second + 1, end); // only one zero-length loop is allowed
    }
  }
  // primary loop moving G-run starting position to the right
  for (s = start; s < end; ++s)
  {
    min_e = s + opts.run_min_len;
    found_any = false;

    // secondary loop moving G-run ending position to the left
    for (e = end; e >= min_e && find_run(s, e, m[i], run_re_c, opts); e--)
    {
      found_any = true;
      // update search bounds
      s = string::const_iterator(m[i].first);
      e = string::const_iterator(m[i].second);
      
      if (opts.fast && i < 3) {
        
        if (m[i].length() == m[i].g_count + 1) {
          // might be a run containing mismatch
          next_tetrad_count = min(tetrad_count, m[i].length());
        } else {
          // perfect or bulged run
          next_tetrad_count = min(tetrad_count, m[i].g_count);
        }
        next_defect_count = defect_count + (m[i].length() != m[i].g_count);
        
        max_score = (next_tetrad_count - 1) * sc.tetrad_bonus
          - next_defect_count * min(sc.bulge_penalty, sc.mismatch_penalty);
        
        if (opts.verbose) {
          print_partial_pqs(m, i, ref);
          Rcout << "next_tetrad_count: " << next_tetrad_count 
                << " next_defect_count: " << next_defect_count
                << " max_scores[0]: " << res.max_scores[m[0].first - ref]
                << " max_score: " << max_score
                << endl;
        }
        if ((max_score < res.max_scores[m[0].first - ref] || max_score < opts.min_score)) {
          // comparison to min_score helps also quite a lot (2-3x speedup for default min_score)
          continue;
        }
      }
      if (i == 0) {
        // enforce G4 total length limit to be relative to the first G-run start
        find_all_runs(
          subject, i+1, e, min(s + opts.max_len, end), m, run_re_c, opts,
          sc, ref, len, pqs_storage, int_cnt, res,
          false, s_time, next_tetrad_count, next_defect_count, fn_call_count
        );
      } else if (i < 3) {
        loop_len = s - m[i-1].second;
        if (loop_len > opts.loop_max_len) {
          return; // skip too long loops
        }
        find_all_runs(
          subject, i+1, e, end, m, run_re_c, opts, sc, ref, len,
          pqs_storage, int_cnt, res,
          (loop_len == 0 ? true : zero_loop), s_time, next_tetrad_count, next_defect_count, fn_call_count
        );
      } else {
        score = 0;
        features_t pqs_features;
        if (opts.use_default_scoring) {
          score = score_pqs(m, pqs_features, sc, opts);
        }
        if ((score || !opts.use_default_scoring) && sc.custom_scoring_fn != NULL) {
          check_custom_scoring_fn(score, m, sc, subject, ref);
        }
        if (score >= opts.min_score) {
          // current PQS satisfied all constraints
          pqs_storage.insert_pqs(score, m[0].first, m[3].second, pqs_features, res);
          
          int pqs_len = m[3].second - m[0].first;
          int offset = m[0].first - ref;
          
          for (int k = 0; k < pqs_len; ++k) {
            res.max_scores[offset + k] = max(res.max_scores[offset + k], score);
            ++res.density[offset + k];
          }
          if (opts.verbose)
            print_pqs(m, score, ref);
        }
      }
    }
    if (!found_any) {
      break;
    }
  }
}


/**
 * Split sequence into individual chunks with sufficient overlap
 * 
 * @param seq DNA sequence
 * @param opts Algorithm options
 */
vector<seq_chunk_t> split_seq_to_chunks(
    const string &seq,
    const opts_t &opts) 
{
  vector<seq_chunk_t> chunk_list;
  
  const size_t seq_len = seq.length();
  const size_t chunk_count = seq_len / opts.chunk_size;
  seq_chunk_t chunk;
  
  for (size_t i = 0; i < chunk_count; ++i) {
    if (i == 0) {
      chunk.s = seq.begin();
    } else {
      chunk.s = seq.begin() + i * opts.chunk_size - 2*opts.max_len;
    }
    if (i == chunk_count - 1) {
      // extend the last chunk to the end of the sequence
      chunk.e = seq.end();
    } else {
      chunk.e = seq.begin() + (i+1) * opts.chunk_size + 2*opts.max_len;
    }
    chunk_list.push_back(chunk);
  }
  if (chunk_count == 0) {
    chunk.s = seq.begin();
    chunk.e = seq.end();
    chunk_list.push_back(chunk);
  }
  return chunk_list;
}


void merge_res_items(results &a, results &b, results &merged) {
  revised_non_overlapping_storage pqs_storage(merged.ref);
  
  // sort results to have them in sequence order
  a.sort_items();
  b.sort_items();
  
  merged.items.clear(); // clear only result items (not vectors)
  vector<results::item_t>::const_iterator left_it = a.items.begin();
  vector<results::item_t>::const_iterator right_it = b.items.begin();
  
  while (left_it != a.items.end() || right_it != b.items.end()) {
    if ((left_it != a.items.end() && right_it == b.items.end()) ||
        (left_it != a.items.end() && right_it != b.items.end() && left_it->start < right_it->start)) {
      pqs_storage.insert_pqs_item(*left_it, merged);
      ++left_it;
    } else {
      pqs_storage.insert_pqs_item(*right_it, merged);
      ++right_it;
    }
  }
  pqs_storage.export_pqs(merged);
}


void find_all_overscored_pqs(
    SEXP subject,
    const string::const_iterator seq_begin,
    const string::const_iterator seq_end,
    const boost::regex &run_re_c,
    const scoring &sc,
    const opts_t &opts,
    results &res,
    int &fn_call_count)
{
  results self_res(seq_end - seq_begin, opts.min_score, seq_begin);
  results a_res(seq_end - seq_begin, opts.min_score, seq_begin);
  results b_res(seq_end - seq_begin, opts.min_score, seq_begin);
  results c_res(seq_end - seq_begin, opts.min_score, seq_begin);
  
  find_overscored_pqs<Overscored::SELF>(
    subject, seq_begin, seq_end, run_re_c, sc, opts, res, self_res, fn_call_count
  );
  // self_res = self overscored
  
  find_overscored_pqs<Overscored::LEFT_BOUND>(
    subject, seq_begin, seq_end, run_re_c, sc, opts, self_res, a_res, fn_call_count
  );
  find_overscored_pqs<Overscored::LEFT_UNBOUND>(
    subject, seq_begin, seq_end, run_re_c, sc, opts, self_res, b_res, fn_call_count
  );
  // a_res = left bound
  // b_res = left unbound
  merge_res_items(a_res, b_res, c_res);
  // c_res = left (un)bound
  
  find_overscored_pqs<Overscored::RIGHT_BOUND>(
    subject, seq_begin, seq_end, run_re_c, sc, opts, self_res, a_res, fn_call_count
  );
  merge_res_items(c_res, a_res, b_res);
  // b_res = left (un)bound + right bound
  
  find_overscored_pqs<Overscored::RIGHT_UNBOUND>(
    subject, seq_begin, seq_end, run_re_c, sc, opts, self_res, a_res, fn_call_count
  );
  merge_res_items(b_res, a_res, c_res);
  // c_res = left (un)bound + right (un)bound
  
  find_overscored_pqs<Overscored::NEIGHBOURING>(
    subject, seq_begin, seq_end, run_re_c, sc, opts, c_res, a_res, fn_call_count
  );
  merge_res_items(c_res, a_res, b_res);
  // b_res = left (un)bound + right (un)bound + neighbouring
  
  find_overscored_pqs<Overscored::NEIGHBOURING_SELF>(
    subject, seq_begin, seq_end, run_re_c, sc, opts, b_res, a_res, fn_call_count
  );
  merge_res_items(b_res, a_res, res);
  // res = left (un)bound + right (un)bound + neighbouring + neighbouring self
}


/**
 * Perform quadruplex search on given DNA sequence.
 *
 * @param subject DNAString object
 * @param seq_begin Beginning of DNA sequence
 * @param seq_end End of DNA sequence
 * @param run_re_c Run regular expression
 * @param sc Scoring options
 * @param opts Algorihtm options
 * @param res Results object
 */
void find_pqs(
    SEXP subject,
    const string::const_iterator seq_begin,
    const string::const_iterator seq_end,
    const boost::regex &run_re_c,
    const scoring &sc,
    const opts_t &opts,
    results &res)
{
  run_match m[RUN_CNT];
  overlapping_storage ov_storage(seq_begin);
  revised_non_overlapping_storage nov_storage(seq_begin);
  storage &pqs_storage = select_pqs_storage(opts.overlapping, ov_storage, nov_storage);
  
  int fn_call_count = 0;
  int int_cnt = 0;
  
  // Global sequence length is the only limit for the first G-run
  find_all_runs(
    subject, 0, seq_begin, seq_end, m, run_re_c, opts, sc,
    seq_begin, seq_end - seq_begin, pqs_storage, int_cnt,
    res, false, chrono::system_clock::now(), INT_MAX, 0, fn_call_count
  );
  pqs_storage.export_pqs(res);
  
  if (opts.fast && !res.items.empty()) {
    size_t prev_count = 0;
    size_t iter_count = 0;
    
    while (prev_count != res.items.size()) {
      ++iter_count;
      prev_count = res.items.size();
      find_all_overscored_pqs(subject, seq_begin, seq_end, run_re_c, sc, opts, res, fn_call_count);
    }
    Rcout << "find_all_overscored_pqs iteration: " << iter_count << endl;
  }
  
  Rcout << "fn_call_count " << fn_call_count << endl;
}


/**
 * Single thread payload
 * 
 * @param tid
 * @param num_threads
 * @param chunk_list
 * @param res_list
 * @param subject
 * @param run_re_c
 * @param sc
 * @param opts
 */
void find_pqs_thread(
    int tid,
    int num_threads,
    vector<seq_chunk_t> &chunk_list,
    vector<results> &res_list,
    SEXP subject,
    const boost::regex &run_re_c,
    const scoring &sc,
    const opts_t &opts)
{
  for (size_t i = tid; i < chunk_list.size(); i += num_threads) {
    find_pqs(
      subject,
      chunk_list[i].s,
      chunk_list[i].e,
      run_re_c,
      sc,
      opts,
      res_list[i]);
  }
}


/**
 * Merge results from different threads
 * TODO: This is probably wrong, because PQS does not go to storage in sequence order
 * 
 * @param seq_begin
 * @param res_list
 * @param res
 * @param opts
 */
void merge_results(
    const string::const_iterator seq_begin,
    vector<results> &res_list,
    results &res,
    const opts_t &opts)
{
  overlapping_storage ov_storage(seq_begin);
  revised_non_overlapping_storage nov_storage(seq_begin);
  storage &pqs_storage = select_pqs_storage(opts.overlapping, ov_storage, nov_storage);
  
  for (size_t i = 0; i < res_list.size(); ++i) {
    // sort partial results to have them in sequence order
    res_list[i].sort_items();
    
    size_t offset = res_list[i].ref - seq_begin;
    
    for (size_t k = 0; k < res_list[i].seq_len; ++k) {
      res.max_scores[offset + k] = max(res.max_scores[offset + k], res_list[i].max_scores[k]);
      res.density[offset + k] = max(res.density[offset + k], res_list[i].density[k]);
    }
    for (size_t k = 0; k < res_list[i].items.size(); ++k) {
      pqs_storage.insert_pqs_item(res_list[i].items[k], res);
    }
  }
  pqs_storage.export_pqs(res);
}


/**
 * Possible split the PQS search between mutliple threads
 * 
 * @param subject
 * @param seq
 * @param run_re_c
 * @param sc
 * @param opts
 * @param res
 */
void find_pqs_parallel(
    SEXP subject,
    const string &seq,
    const boost::regex &run_re_c,
    const scoring &sc,
    const opts_t &opts,
    results &res)
{
  if (opts.threads == 1) {
    find_pqs(subject, seq.begin(), seq.end(), run_re_c, sc, opts, res);
  } else {
    vector<seq_chunk_t> chunk_list = split_seq_to_chunks(seq, opts);
    
    size_t num_threads = min(chunk_list.size(), opts.threads);
    boost::thread *tt;
    
    boost::thread::attributes attrs;
    
    // initialize result objects
    vector<results> res_list(chunk_list.size());
    for (size_t i = 0; i < chunk_list.size(); ++i) {
      res_list[i].init(chunk_list[i].e - chunk_list[i].s, opts.min_score, chunk_list[i].s);
    }
    if (num_threads > 1) {
      // run additional threads
      tt = new boost::thread[num_threads - 1];
      for (size_t tid = 1; tid < num_threads; ++tid) {
        tt[tid - 1] = boost::thread(
          attrs, boost::bind(find_pqs_thread, tid, num_threads, std::ref(chunk_list), std::ref(res_list),
          subject, std::ref(run_re_c), std::ref(sc), std::ref(opts))
        );
      }
    }
    // run main thread calculation
    find_pqs_thread(0, num_threads, chunk_list, res_list, subject, run_re_c, sc, opts);
    
    if (num_threads > 1) {
      // wait for additional threads
      for (size_t tid = 1; tid < num_threads; ++tid) {
        tt[tid - 1].join();
      }
      delete [] tt;
    }
    merge_results(seq.begin(), res_list, res, opts);
  }
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
//' @param threads Number of threads that can be created by pqsfinder to
//'   speed up the computation.
//' @param fast Enable fast searching. This has some impact on maxScores and
//'   density vectors.
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
    bool fast = true,
    int threads = 1,
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
  
  opts_t opts;
  opts.use_re = false;
  opts.debug = false;
  opts.verbose = verbose;
  opts.use_default_scoring = use_default_scoring;
  opts.fast = fast;
  opts.overlapping = overlapping;
  opts.max_len = max_len;
  opts.min_score = min_score;
  opts.loop_max_len = loop_max_len;
  opts.loop_min_len = loop_min_len;
  opts.run_max_len = run_max_len;
  opts.run_min_len = run_min_len;
  opts.threads = threads;
  opts.chunk_size = 500 * max_len;
  
  if (threads > 1) {
    // disable verbose mode
    opts.verbose = false;
  }
  if (run_re != "G{1,10}.{0,9}G{1,10}") {
    // User specified its own regexp, force to use regexp engine
    opts.use_re = true;
  }
  if (opts.use_re) {
    opts.check_int_period = 1e6;
  } else {
    opts.check_int_period = 1e6;
  }
  if (opts.overlapping) {
    // cannot use optimization when searching for overlapping G4s
    opts.fast = false;
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
  sc.init_loop_penalties(opts.loop_max_len);
  sc.init_bulge_penalties(opts.run_max_len);

  string seq = as<string>(as_character(subject));
  Function reverseComplement("reverseComplement");
  SEXP subject_rc = reverseComplement(subject);
  string seq_rc = as<string>(as_character(subject_rc));

  results res_sense(seq.length(), opts.min_score, seq.begin());
  results res_antisense(seq_rc.length(), opts.min_score, seq_rc.begin());
  boost::regex run_re_c(run_re);

  if (opts.debug) {
    Rcout << "G-run regexp: " << run_re << endl;
    Rcout << "Use regexp engine: " << opts.use_re << endl;
    Rcout << "Input sequence length: " << seq.length() << endl;
    Rcout << "Use user fn: " << (custom_scoring_fn != R_NilValue) << endl;
  }

  if (custom_scoring_fn != R_NilValue) {
    sc.set_custom_scoring_fn(custom_scoring_fn);
    if (opts.debug) {
      Rcout << "User function: " << endl;
      Function show("show");
      show(custom_scoring_fn);
    }
  }

  #ifdef GPERF_ENABLED
    ProfilerStart("profiling.log");
  #endif

  if (strand == "+" || strand == "*") {
    Rcout << "Searching on sense strand..." << endl;
    find_pqs_parallel(subject, seq, run_re_c, sc, opts, res_sense);
    Rcout << "Search status: finished              " << endl;
  }
  if (strand == "-" || strand == "*") {
    Rcout << "Searching on antisense strand..." << endl;
    find_pqs_parallel(subject_rc, seq_rc, run_re_c, sc, opts, res_antisense);
    Rcout << "Search status: finished              " << endl;
  }

  #ifdef GPERF_ENABLED
    ProfilerStop();
  #endif
  
  size_t res_items_size = res_sense.items.size() + res_antisense.items.size();

  IntegerVector res_start(res_items_size);
  IntegerVector res_width(res_items_size);
  IntegerVector res_score(res_items_size);
  CharacterVector res_strand(res_items_size);
  IntegerVector res_nt(res_items_size);
  IntegerVector res_nb(res_items_size);
  IntegerVector res_nm(res_items_size);
  IntegerVector res_rl1(res_items_size);
  IntegerVector res_rl2(res_items_size);
  IntegerVector res_rl3(res_items_size);
  IntegerVector res_ll1(res_items_size);
  IntegerVector res_ll2(res_items_size);
  IntegerVector res_ll3(res_items_size);
  
  size_t i, k;
  
  for (i = 0; i < res_sense.items.size(); ++i) {
    // s - ref + 1; R indexing starts at 1
    res_start[i] = res_sense.items[i].start - seq.begin() + 1;
    res_width[i] = res_sense.items[i].len;
    res_score[i] = res_sense.items[i].score;
    res_strand[i] = "+";
    res_nt[i] = res_sense.items[i].nt;
    res_nb[i] = res_sense.items[i].nb;
    res_nm[i] = res_sense.items[i].nm;
    res_rl1[i] = res_sense.items[i].rl1;
    res_rl2[i] = res_sense.items[i].rl2;
    res_rl3[i] = res_sense.items[i].rl3;
    res_ll1[i] = res_sense.items[i].ll1;
    res_ll2[i] = res_sense.items[i].ll2;
    res_ll3[i] = res_sense.items[i].ll3;
  }
  
  for (k = 0; k < res_antisense.items.size(); ++k) {
    // seq_len - (e - ref) + 1; R indexing starts at 1
    res_start[i] = seq_rc.length() -
      (res_antisense.items[k].start +
      res_antisense.items[k].len -
      seq_rc.begin()) + 1;
    res_width[i] = res_antisense.items[k].len;
    res_score[i] = res_antisense.items[k].score;
    res_strand[i] = "-";
    res_nt[i] = res_antisense.items[k].nt;
    res_nb[i] = res_antisense.items[k].nb;
    res_nm[i] = res_antisense.items[k].nm;
    res_rl1[i] = res_antisense.items[k].rl1;
    res_rl2[i] = res_antisense.items[k].rl2;
    res_rl3[i] = res_antisense.items[k].rl3;
    res_ll1[i] = res_antisense.items[k].ll1;
    res_ll2[i] = res_antisense.items[k].ll2;
    res_ll3[i] = res_antisense.items[k].ll3;
    ++i;
  }

  IntegerVector res_density(seq.length());
  IntegerVector res_max_scores(seq.length());
  
  for (size_t i = 0; i < seq.length(); ++i) {
    size_t antisense_i = seq.length() - i - 1;
    res_density[i] = res_sense.density[i] + res_antisense.density[antisense_i];
    res_max_scores[i] = max(res_sense.max_scores[i], res_antisense.max_scores[antisense_i]);
  }
  Function pqsviews("PQSViews");
  
  return pqsviews(
    subject, res_start, res_width, res_strand, res_score,
    res_density, res_max_scores,
    res_nt, res_nb, res_nm,
    res_rl1, res_rl2, res_rl3,
    res_ll1, res_ll2, res_ll3);
}
