/**
 * Run match structure
 *
 * Author: Jiri Hon <jiri.hon@gmail.com>
 * Date: 2018/02/11
 * Package: pqsfinder
 */

#ifndef RUN_MATCH_HEADER
#define RUN_MATCH_HEADER

using namespace std;

static const int RUN_CNT = 4;


class run_match {
public:
  string::const_iterator first;
  string::const_iterator second;
  int g_count;
  int length() const {
    return second - first;
  };
};


#endif // RUN_MATCH_HEADER
