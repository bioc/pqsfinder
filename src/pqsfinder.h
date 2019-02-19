/**
 * pqsfinder header file
 *
 * Author: Jiri Hon <jiri.hon@gmail.com>
 * Date: 2018/02/11
 * Package: pqsfinder
 */

#ifndef PQSFINDER_HEADER
#define PQSFINDER_HEADER

#include "boost/regex.hpp"
#include "opts.h"
#include "results.h"
#include "run_match.h"
#include "scoring.h"
#include "storage.h"

using namespace std;


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
    int &pqs_cnt,
    results &res,
    bool zero_loop,
    const chrono::system_clock::time_point s_time,
    int tetrad_count,
    int defect_count,
    // int loop_sum,
    int &fn_call_count,
    bool show_progress);


#endif // PQSFINDER_HEADER
