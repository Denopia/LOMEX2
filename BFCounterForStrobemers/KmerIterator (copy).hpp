#ifndef BFG_KMER_ITERATOR_HPP
#define BFG_KMER_ITERATOR_HPP

#include <iterator>
#include <iostream>
#include <vector>
#include <deque>
#include <map>
#include <unordered_map>
#include <tuple>
#include "Kmer.hpp"


/* Short description: 
 *  - Easily iterate through kmers in a read
 *  - If the read contains any N, then the N is skipped and checked whether
 *    there is a kmer to the right of the N
 *  - iter->first gives the kmer, iter->second gives the position within the reads
 * */
class KmerIterator : public std::iterator<std::input_iterator_tag, std::pair<Kmer, int>, int> {
public:
  KmerIterator() : s_(NULL), p_(), invalid_(true), total_length_(0), care_length_(0), care_status_(NULL) {}
  KmerIterator(const char* s) : s_(s), p_(), invalid_(false), total_length_(0), care_length_(0), care_status_(NULL) {find_next(-1,-1,false);}

  KmerIterator(const char* s, int tl, int cl, std::vector<bool> * cs) : s_(s), p_(), invalid_(false), total_length_(tl), care_length_(cl), care_status_(cs) 
  {
    p_.second = -1;
    while(!find_next_II(p_.second, p_.second+total_length_, false)){};
  }

  KmerIterator(const char* s, int mw, int mn, int ml) : s_(s), p_(), invalid_(false), min_window_(mw), min_num_(mn), min_len_(ml) 
  {
    p_.second = 0;
    previous_minimizer_kmer_ = "IAMEMPTY";
    found_minimizer_kmers_.clear();
    //while(!find_next_III(p_.second)){};
    //std::cout << "HOLA\n"; 
    //find_minimizers();
    //find_next_III(0);
  }

  KmerIterator(const KmerIterator& o) : s_(o.s_), p_(o.p_), invalid_(o.invalid_), total_length_(o.total_length_), care_length_(o.care_length_), care_status_(o.care_status_) {}


  //KmerIterator() : s_(NULL), p_(), invalid_(true), total_length_(0), care_length_(0), care_status_(NULL) {std::cout << "LOOKS BAD" << std::endl;}
  //KmerIterator(const char* s) : s_(s), p_(), invalid_(false), total_length_(0), care_length_(0), care_status_(NULL) { std::cout << "LOOKS NO GOOD" << std::endl;find_next(-1,-1,false);}
  //KmerIterator(const char* s, int tl, int cl, std::vector<bool> * cs) : s_(s), p_(), invalid_(false), total_length_(tl), care_length_(cl), care_status_(cs) { std::cout << "LOOKS GOOD" << std::endl;find_next_II(-1,-1,false);}
  //KmerIterator(const KmerIterator& o) : s_(o.s_), p_(o.p_), invalid_(o.invalid_), total_length_(o.total_length_), care_length_(o.care_length_), care_status_(o.care_status_) {std::cout << "LOOKS GOOD" << std::endl;}


  KmerIterator& operator++();
  KmerIterator operator++(int);
  void raise(Kmer &km, Kmer &rep);

  bool operator==(const KmerIterator& o);
  bool operator!=(const KmerIterator& o) { return !this->operator==(o);}

  std::pair<Kmer, int>& operator*();
  std::pair<Kmer, int>* operator->();

  //std::tuple<Kmer, bool> extract_spaced_kmer();
  bool extract_spaced_kmer();
  void find_minimizers();
  void find_minimizers_faster();
  void find_minimizers_faster_super();
  void minimizer_check(int index);

  bool imfine;

private:
  void find_next(size_t i, size_t j, bool last_valid);
  bool find_next_II(size_t i, size_t j, bool last_valid);
  bool find_next_III();
  uint64_t minimizer_hash(uint64_t key);
  uint64_t string2int64(int start, int strand);
  char reverse_complement_character(char c);


  // New stuff for the non-minimizer version
  std::vector<bool> * care_status_;
  unsigned int total_length_;
  unsigned int care_length_;
  unsigned int window_length_;

  // These are for the minimizer version

  std::map<int, std::vector<std::tuple<uint64_t,int,int> > > minimizers_; // map with minimizer window position as key, value is tuple with [hash value, position, strand(0=forward, 1=reverse)]
  std::map<int, std::vector<std::tuple<uint64_t,int,int> > > minimizers_faster_; // map with minimizer window position as key, value is tuple with [hash value, position, strand(0=forward, 1=reverse)]

  std::map<int, std::vector<std::tuple<uint64_t,int,int> > > minimizers_faster_super_; // map with minimizer window position as key, value is tuple with [hash value, position, strand(0=forward, 1=reverse)]
  std::map<int, std::vector<std::tuple<uint64_t,int,int> > > minimizers_faster_super_reverse; // map with minimizer window position as key, value is tuple with [hash value, position, strand(0=forward, 1=reverse)]
  
  unsigned int read_len_; // read length
  unsigned int min_len_; // minimizer length
  unsigned int min_window_; // number of minimizers in a window  
  unsigned int min_num_; // number of minimizers used for one spaced k-mer
  std::vector<std::tuple<int,int> > minimizer_kmer_fixed_streak_starts_;
  

  std::string previous_minimizer_kmer_;
  std::unordered_map<std::string, int> found_minimizer_kmers_;
  
  //std::string kmer_window_; // unnecessary

  const char *s_;
  //Kmer spaced_kmer;
  //std::pair<Kmer, int> P_; // For the whole window
  std::deque<char> window_;
  std::deque<bool> window_status_; // For validity indicators 
  //std::pair<Kmer, int> v_; 
  std::pair<Kmer, int> p_; // For the spaced k-mers
  bool invalid_;
};

#endif // BFG_KMER_ITERATOR_HPP
