/**
 * Storage class aggregating locally-best PQS.
 *
 * Author: Jiri Hon <jiri.hon@gmail.com>
 * Date: 2016/03/30
 * Package: pqsfinder
 */

#ifndef STORAGE_HEADER
#define STORAGE_HEADER

#include <Rcpp.h>
#include "features.h"
#include "opts.h"
#include "results.h"

using namespace Rcpp;
using namespace std;


class storage {
public:
  virtual ~storage() {};
  virtual void insert_pqs(
      int score, string::const_iterator s, string::const_iterator e, features_t &f,
      results &res) = 0;
  virtual void insert_pqs_item(const results::item_t &item, results &res) {
    features_t f;
    f.nt = item.nt;
    f.nb = item.nb;
    f.nm = item.nm;
    f.rl1 = item.rl1;
    f.rl2 = item.rl2;
    f.rl3 = item.rl3;
    f.ll1 = item.ll1;
    f.ll2 = item.ll2;
    f.ll3 = item.ll3;
    
    this->insert_pqs(
      item.score,
      item.start,
      item.start + item.len,
      f, res);
  }
  virtual void export_pqs(results &res) = 0;
};

class overlapping_storage : public storage {
private:
  typedef struct {
    int score;
    features_t f;
  } map_value_t;
  typedef map< string::const_iterator, map_value_t > map_t;
  map_t pqs_map;
  string::const_iterator pqs_start;
  
public:
  overlapping_storage(string::const_iterator pqs_start) : pqs_start(pqs_start) {}
  
  virtual void insert_pqs(
      int score, string::const_iterator s, string::const_iterator e, features_t &f,
      results &res)
  {
    if (this->pqs_start < s) {
      this->export_pqs(res); // this clears the pqs_map
      this->pqs_start = s;
    }
    map_t::iterator it = this->pqs_map.find(e);
    if (it != this->pqs_map.end()) {
      if (score > it->second.score) {
        it->second.score = score;
        it->second.f = f;
      } else {
        return;
      }
    } else {
      map_value_t value = {score, f};
      this->pqs_map.insert(make_pair(e, value));
    }
  }
  virtual void export_pqs(
      results &res)
  {
    for (map_t::iterator it = this->pqs_map.begin(); it != this->pqs_map.end(); ++it) {
      res.save_pqs(it->second.score, this->pqs_start, it->first, it->second.f);
    }
    this->pqs_map.clear();
  }
};

class revised_non_overlapping_storage: public storage {
private:
  class range {
  public:
    string::const_iterator s;
    string::const_iterator e;
    features_t f;
    range(string::const_iterator &s, string::const_iterator &e, features_t &f) :
      s(s), e(e), f(f) {};
    range() {};
  };
  typedef map< int, list<range> > storage_t; // map is always sorted by keys
  storage_t st;
  string::const_iterator last_e;
  
public:
  revised_non_overlapping_storage(string::const_iterator start) : last_e(start) {}
  
  virtual void insert_pqs(
      int score, string::const_iterator s, string::const_iterator e, features_t &f,
      results &res)
  {
    if (s >= this->last_e && !this->st.empty())
    {// export PQS because no further overlapping pqs can be found
      this->export_pqs(res);
    }
    if (e > this->last_e) {
      this->last_e = e;
    }
    storage_t::iterator it = this->st.find(score);
    if (it != this->st.end()) {
      list<range> &list = it->second;
      if (list.empty()) {
        throw runtime_error("Inconsistent state of non-overlapping storage.");
      }
      // Insert new pqs
      list.push_back(range(s, e, f));
    }
    else {
      this->st.insert(storage_t::value_type(score, list<range>(1, range(s, e, f))));
    }
  }
  virtual void export_pqs(
      results &res)
  {
    range best_pqs;
    storage_t::iterator it, r_it, temp;
    list<range>::iterator prev, curr;
    
    while (!this->st.empty()) {
      it = --this->st.end(); // decrement to point on last list
      // [it] points on the highest scoring list of ranges
      
      // resolve overlaps between equal-scoring PQS
      prev = it->second.begin();
      curr = next(prev);
      
      while (curr != it->second.end()) {
        // if (prev->s <= curr->s && curr->e <= prev->e) {
        //   it->second.erase(prev);
        //   prev = curr;
        //   ++curr;
        // } else if (prev->e > curr->s) {
        //   it->second.erase(curr);
        //   curr = next(prev);
        // } else {
        //   prev = curr;
        //   ++curr;
        // }
        if (curr->s < prev->e) {
          if ((curr->e - curr->s) < (prev->e - prev->s)) {
            // curr is shorter
            it->second.erase(prev);
            prev = curr;
            ++curr;
          } else {
            // prev is shorter or equal, so leave it
            it->second.erase(curr);
            curr = next(prev);
          }
        } else {
          prev = curr;
          ++curr;
        }
      }
      if (it->second.empty()) {
        throw runtime_error("Inconsistent state of non-overlapping PQS list.");
      }
      while (!it->second.empty()) {
        best_pqs = it->second.front();
        res.save_pqs(it->first, best_pqs.s, best_pqs.e, best_pqs.f);
        it->second.pop_front();
        
        if (it != this->st.begin()) {
          r_it = std::prev(it); // set score level to the next lower level
          
          while (true) {
            // remove all overlapping PQS with lower score
            list<range> &l = r_it->second;
            list<range>::iterator l_it = l.begin(), l_temp;
            
            while (l_it != l.end()) {
              if ((l_it->s <= best_pqs.s && best_pqs.s < l_it->e) ||
                  (best_pqs.s <= l_it->s && l_it->s < best_pqs.e)) {
                l_temp = l_it; // erase operation invalidates iterator
                ++l_it;
                l.erase(l_temp);
              } else {
                ++l_it;
              }
            }
            if (r_it == this->st.begin()) {
              // the end of iteration
              if (l.empty()) {
                this->st.erase(r_it); // erase empty score level
              }
              break;
            } else if (l.empty()) {
              // delete empty score level from storage and move on lower score level
              temp = r_it; // erase operation invalidates iterator
              --r_it;
              this->st.erase(temp);
            }
            else {
              // just move on lower score level
              --r_it;
            }
          }
        }
      }
      this->st.erase(it); // erase empty score level
    }
  }
  inline void print() {
    for (storage_t::const_iterator it = this->st.begin(); it != this->st.end(); ++it) {
      Rcout << it->first << ":";
      for (list<range>::const_iterator lit = it->second.begin(); lit != it->second.end(); ++lit) {
        Rcout << " " << string(lit->s, lit->e);
      }
      Rcout << endl;
    }
  }
};


class fast_non_overlapping_storage: public storage {
private:
  string::const_iterator best_s, best_e;
  features_t best_f;
  int best_score = 0;
  
public:
  fast_non_overlapping_storage(string::const_iterator start) : best_e(start) {}
  
  virtual void insert_pqs(
      int score, string::const_iterator s, string::const_iterator e, features_t &f,
      results &res)
  {
    if (s >= this->best_e && this->best_score > 0)
    {// export PQS because no further overlapping pqs can be found
      this->export_pqs(res);
    }
    if (score > this->best_score ||
        (score == this->best_score && (e - s) < (this->best_e - this->best_s))) {
      
      int pqs_len = e - s;
      int offset = s - res.ref;
      
      for (int k = 0; k < pqs_len; ++k) {
        res.max_scores[offset + k] = max(res.max_scores[offset + k], score);
      }
      
      if (e < this->best_e) {
        // reset max scores on the right side of the new best
        int len = this->best_e - e;
        int offset = e - res.ref;
        
        for (int k = 0; k < len; ++k) {
          res.max_scores[offset + k] = 0;
        }
      }
      this->best_score = score;
      this->best_s = s;
      this->best_e = e;
      this->best_f = f;
      // Rcout << "new best: " << s - res.ref + 1 << "-" << e - res.ref << " " << score << endl;
    }
  }
  virtual void export_pqs(
      results &res)
  {
    if (this->best_score > 0) {
      res.save_pqs(this->best_score, this->best_s, this->best_e, this->best_f);
      this->best_score = 0;
    }
  }
};


/**
 * Select between overlapping and non-overlapping storage.
 * 
 * @param opts Algorithm options.
 * @param ov Overlapping storage.
 * @param nov Non-overlapping storage.
 * @return Reference to storage interface.
 */
storage &select_pqs_storage(
    const opts_t &opts,
    overlapping_storage &ov,
    revised_non_overlapping_storage &nov,
    fast_non_overlapping_storage &gnov)
{
  if (opts.overlapping) {
    return ov;
  } else if (opts.fast) {
    return gnov;
  } else {
    return nov;
  }
}

#endif // STORAGE_HEADER
