/**
 * Storage class aggregating locally-best PQS.
 *
 * Author: Jiri Hon <jiri.hon@gmail.com>
 * Date: 2016/03/30
 * Package: pqsfinder
 */

#ifndef PQS_STORAGE_HEADER
#define PQS_STORAGE_HEADER

#include <Rcpp.h>
#include "results.h"

using namespace Rcpp;
using namespace std;

class pqs_storage {
public:
  virtual ~pqs_storage() {};
  virtual void insert_pqs(
      int score, string::const_iterator s, string::const_iterator e,
      results &res, const string::const_iterator &ref, const string &strand) = 0;
  virtual void export_pqs(
      results &res, const string::const_iterator &ref,
      const string &strand) = 0;
};

class pqs_storage_overlapping : public pqs_storage {
private:
  typedef map< string::const_iterator, int > map_t;
  map_t pqs_map;
  string::const_iterator pqs_start;
  
public:
  pqs_storage_overlapping(string::const_iterator pqs_start) : pqs_start(pqs_start) {}
  
  virtual void insert_pqs(
      int score, string::const_iterator s, string::const_iterator e,
      results &res, const string::const_iterator &ref, const string &strand)
  {
    if (this->pqs_start < s) {
      this->export_pqs(res, ref, strand);
      this->pqs_start = s;
    }
    map_t::iterator it = pqs_map.find(e);
    if (it != pqs_map.end()) {
      if (score > it->second) {
        it->second = score;
      } else {
        return;
      }
    } else {
      pqs_map.insert(make_pair(e, score));
    }
  }
  virtual void export_pqs(
      results &res, const string::const_iterator &ref,
      const string &strand)
  {
    for (map_t::iterator it = this->pqs_map.begin(); it != this->pqs_map.end(); ++it) {
      res.save_pqs(it->second, this->pqs_start, it->first, ref, strand);
    }
    pqs_map.clear();
  }
};

class pqs_storage_non_overlapping: public pqs_storage {
private:
  class range {
  public:
    string::const_iterator s;
    string::const_iterator e;
    range(string::const_iterator s, string::const_iterator e) :
      s(s), e(e) {};
    range() {};
  };
  typedef map< int, list<range> > storage_t;
  storage_t st;
  
  typedef struct pqs {
    int score;
    string::const_iterator s;
    string::const_iterator e;
  } pqs_t;
  pqs_t best;

public:
  pqs_storage_non_overlapping() {
    this->best.score = 0;
  }
  virtual void insert_pqs(
      int score, string::const_iterator s, string::const_iterator e,
      results &res, const string::const_iterator &ref, const string &strand)
  {
    if (this->best.score && s >= this->best.e)
    {// export PQS because no further overlapping pqs can be found
      this->export_pqs(res, ref, strand);
    }
    if (score > this->best.score ||
        (score == this->best.score &&
         this->best.s <= s && e <= this->best.e)) {
      this->best.score = score;
      this->best.s = s;
      this->best.e = e;
    }
    storage_t::iterator it = st.find(score);
    if (it != st.end()) {
      list<range> &list = it->second;
      if (list.empty()) {
        list.push_back(range(s, e));
      }
      else {
        range &last = list.back();
        if (last.s <= s && e <= last.e) {
          // Replace existing by shorter one
          last.s = s;
          last.e = e;
        }
        else if (last.e <= s) {
          // Insert new non-overlapping pqs
          list.push_back(range(s, e));
        }
      }
    }
    else {
      st.insert(storage_t::value_type(score, list<range>(1, range(s, e))));
    }
  }
  virtual void export_pqs(
      results &res, const string::const_iterator &ref,
      const string &strand)
  {
    this->best.score = 0; // reset

    range best_pqs;
    storage_t::iterator it;
    storage_t::iterator temp;

    while (!st.empty()) {
      it = --st.end(); // decrement to point on last list
      best_pqs = it->second.back();

      res.save_pqs(it->first, best_pqs.s, best_pqs.e, ref, strand);

      while (true) {
        // remove all overlapping PQS with lower score
        // Rcout << "Removing from list " << it->first << endl;
        list<range> &list = it->second;
        while (!list.empty() &&
               ((list.back().s <= best_pqs.s && best_pqs.s < list.back().e) ||
               (best_pqs.s <= list.back().s && list.back().s < best_pqs.e)))  {
          list.pop_back();
        }
        if (it == st.begin()) {
          // the end of iteration
          if (list.empty())
            st.erase(it);
          break;
        }
        else if (list.empty()) {
          // delete empty list from storage
          temp = it; // erase operation invalidates iterator
          --it;
          st.erase(temp);
        }
        else {
          --it;
        }
      }
    }
  }
  inline void print() {
    for (storage_t::const_iterator it = st.begin(); it != st.end(); ++it) {
      Rcout << it->first << ":";
      for (list<range>::const_iterator lit = it->second.begin(); lit != it->second.end(); ++lit) {
        Rcout << " " << string(lit->s, lit->e);
      }
      Rcout << endl;
    }
  }
};

#endif // PQS_STORAGE_HEADER
