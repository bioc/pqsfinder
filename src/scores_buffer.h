/**
 * Max scores buffer
 *
 * Author: Jiri Hon <jiri.hon@gmail.com>
 * Date: 2018/02/11
 * Package: pqsfinder
 */

#ifndef SCORES_BUFFER_HEADER
#define SCORES_BUFFER_HEADER

using namespace std;

class scores_buffer {
private:
  const int max_len;
  const string::const_iterator ref;
  int *buffer = NULL;
  string::const_iterator bstart;
  int bstart_offset;
  
public:
  scores_buffer(const int max_len, const string::const_iterator ref):
    max_len(max_len), ref(ref), bstart(ref), bstart_offset(0)  {
    
    this->buffer = new int[max_len];
    for (int i = 0; i < max_len; ++i) {
      this->buffer[i] = 0;
    }
  }
  ~scores_buffer() {
    delete [] this->buffer;
  }
  inline void move_and_set(
      const string::const_iterator s, const string::const_iterator e, const int score)
  {
    int items_to_reset;
    int diff = s - this->bstart;
    
    if (diff < 0) {
      // clear whole buffer when moving back;
      items_to_reset = this->max_len;
    }
    items_to_reset = min(diff, this->max_len);
    
    if (items_to_reset > 0) {
      int curr_offset = this->bstart - this->ref;
      
      for (int i = 0; i < items_to_reset; ++i) {
        int k = (curr_offset + i) % this->max_len;
        this->buffer[k] = 0;
      }
    }
    // set buffer start pointer and offset
    this->bstart = s;
    this->bstart_offset = s - this->ref;
    
    // set max scores
    int write_len = e - s;
    for (int i = 0; i < write_len; ++i) {
      int k = (this->bstart_offset + i) % this->max_len;
      this->buffer[k] = max(this->buffer[k], score);
    }
  }
  inline int get(const string::const_iterator s) {
    int diff = s - this->bstart;
    if (diff < this->max_len) {
      int offset = (this->bstart_offset + diff) % this->max_len;
      return this->buffer[offset];
    } else {
      return 0;
    }
  }
  inline void clear_range(const string::const_iterator s, const string::const_iterator e) {
    int reset_len = e - s;
    int reset_offset = s - this->bstart;
    int max_items_to_reset = max(this->max_len - reset_offset, 0);
    int items_to_reset = min(reset_len, max_items_to_reset);
    
    int offset = s - this->ref;
    
    for (int i = 0; i < items_to_reset; ++i) {
      int k = (offset + i) % this->max_len;
      this->buffer[k] = 0;
    }
  }
  inline void clear() {
    for (int i = 0; i < this->max_len; ++i) {
      this->buffer[i] = 0;
    }
    this->bstart = this->ref;
    this->bstart_offset = 0;
  }
};

#endif // SCORES_BUFFER_HEADER
