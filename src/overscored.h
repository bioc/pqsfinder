/**
 * Functions to find overscored PQS when using fast heuristics
 *
 * Author: Jiri Hon <jiri.hon@gmail.com>
 * Date: 2018/02/11
 * Package: pqsfinder
 */

#ifndef OVERSCORED_HEADER
#define OVERSCORED_HEADER

#include <Rcpp.h>
#include <boost/regex.hpp>
#include "opts.h"
#include "pqsfinder.h"
#include "results.h"
#include "run_match.h"
#include "scoring.h"
#include "storage.h"

namespace Overscored {
enum type_t {
  SELF,
  LEFT_BOUND,
  LEFT_UNBOUND,
  RIGHT_BOUND,
  RIGHT_UNBOUND,
  NEIGHBOURING,
  NEIGHBOURING_SELF,
};
}


/**
 * Find PQS that were missed during fast search
 * 
 * @param subject
 * @param seq_begin Beginning of DNA sequence
 * @param seq_end End of DNA sequence
 * @param run_re_c
 * @param sc
 * @param opts
 * @param res
 * @param new_res
 * @param pqs_storage
 * @param fn_call_count
 */
template <Overscored::type_t type>
void find_overscored_pqs(
    SEXP subject,
    const string::const_iterator seq_begin,
    const string::const_iterator seq_end,
    const boost::regex &run_re_c,
    const scoring &sc,
    const opts_t &opts,
    results &res,
    results &new_res,
    int &fn_call_count
)
{
  run_match m[RUN_CNT];
  int pqs_cnt = 0;
  string::const_iterator left_start, left_end, right_start, right_end,
    next_pqs_start, prev_pqs_end, search_start, search_end, left_unbound, right_unbound;
  
  res.sort_items(); // sort results to have them in sequence order
  new_res.items.clear();
  
  for (size_t i = 0; i < res.items.size(); ++i) {
    
    // new_res.clear_max_scores();
    
    revised_non_overlapping_storage pqs_storage(seq_begin);
    
    left_end = res.items[i].start;
    right_start = res.items[i].start + res.items[i].len;
    
    if (i == 0) {
      left_start = max(left_end - opts.max_len, seq_begin);
    } else {
      prev_pqs_end = res.items[i-1].start + res.items[i-1].len;
      left_start = max(left_end - opts.max_len, prev_pqs_end);
      if (type == Overscored::NEIGHBOURING && left_start - prev_pqs_end < opts.max_len) {
        left_start = left_end; // do not search again
      }
    }
    if (i == res.items.size() - 1) {
      right_end = min(right_start + opts.max_len, seq_end);
    } else {
      next_pqs_start = res.items[i+1].start;
      right_end = min(right_start + opts.max_len, next_pqs_start);
      if (type == Overscored::NEIGHBOURING && next_pqs_start - right_end < opts.max_len) {
        right_end = next_pqs_start; // extend search region
      }
    }
    left_unbound = max(left_end - opts.max_len, seq_begin);
    right_unbound = min(right_start + opts.max_len, seq_end);
    
    // clear max scores in the region
    for (int i = left_unbound - seq_begin; i < right_unbound - seq_begin; ++i) {
      new_res.max_scores[i] = 0;
    }
    
    if (type == Overscored::SELF) {
      find_all_runs(
        subject, 0, left_end, right_start, m, run_re_c, opts, sc, 
        seq_begin, seq_end - seq_begin, pqs_storage,
        pqs_cnt, new_res, false, chrono::system_clock::now(),
        INT_MAX, 0, fn_call_count
      );
    } else if (type == Overscored::NEIGHBOURING) {
      if (left_end - left_start > opts.run_min_len * 4) {
        find_all_runs(
          subject, 0, left_start, left_end, m, run_re_c, opts, sc, 
          seq_begin, seq_end - seq_begin, pqs_storage,
          pqs_cnt, new_res, false, chrono::system_clock::now(),
          INT_MAX, 0, fn_call_count
        );
      }
      if (right_end - right_start > opts.run_min_len * 4) {
        find_all_runs(
          subject, 0, right_start, right_end, m, run_re_c, opts, sc, 
          seq_begin, seq_end - seq_begin, pqs_storage,
          pqs_cnt, new_res, false, chrono::system_clock::now(),
          INT_MAX, 0, fn_call_count
        );
      }
    } else if (type == Overscored::LEFT_BOUND) {
      // Rcout << "region " << left_start - seq_begin + 1 << "-" << right_start - seq_begin << endl;
      find_all_runs(
        subject, 0, left_start, right_start, m, run_re_c, opts, sc, 
        seq_begin, seq_end - seq_begin, pqs_storage,
        pqs_cnt, new_res, false, chrono::system_clock::now(),
        INT_MAX, 0, fn_call_count
      );
    } else if (type == Overscored::LEFT_UNBOUND) {
      find_all_runs(
        subject, 0, left_unbound, right_start, m, run_re_c, opts, sc, 
        seq_begin, seq_end - seq_begin, pqs_storage,
        pqs_cnt, new_res, false, chrono::system_clock::now(),
        INT_MAX, 0, fn_call_count
      );
    } else if (type == Overscored::RIGHT_BOUND) {
      find_all_runs(
        subject, 0, left_end, right_end, m, run_re_c, opts, sc, 
        seq_begin, seq_end - seq_begin, pqs_storage,
        pqs_cnt, new_res, false, chrono::system_clock::now(),
        INT_MAX, 0, fn_call_count
      );
    } else if (type == Overscored::RIGHT_UNBOUND) {
      // Rcout << "region " << left_end - seq_begin + 1 << "-" <<  min(right_start + opts.max_len, seq_end) - seq_begin << endl;
      find_all_runs(
        subject, 0, left_end, right_unbound, m, run_re_c, opts, sc, 
        seq_begin, seq_end - seq_begin, pqs_storage,
        pqs_cnt, new_res, false, chrono::system_clock::now(),
        INT_MAX, 0, fn_call_count
      );
    } else if (type == Overscored::NEIGHBOURING_SELF) {
      find_all_runs(
        subject, 0, left_start, right_end, m, run_re_c, opts, sc, 
        seq_begin, seq_end - seq_begin, pqs_storage,
        pqs_cnt, new_res, false, chrono::system_clock::now(),
        INT_MAX, 0, fn_call_count
      );
    }
    pqs_storage.export_pqs(new_res);
    
    // if (type == Overscored::RIGHT_UNBOUND) {
      // new_res.sort_items();
      // new_res.print();
    // }
  }
  // Rcout << "After solve overscored " << type << endl;
  // new_res.sort_items();
  // new_res.print();
}

#endif // OVERSCORED_HEADER
