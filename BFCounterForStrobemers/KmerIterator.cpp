#include <iterator>
#include <utility>
#include <limits>
#include <queue>
#include "Kmer.hpp"
#include "KmerIterator.hpp"


/* Note: That an iter is exhausted means that (iter._invalid == true) */

// use:  ++iter;
// pre:  
// post: *iter is now exhausted
//       OR *iter is the next valid pair of kmer and location
KmerIterator& KmerIterator::operator++() {

  //std::cout << "ITERATOR ++ V1 OPERATION STARTS\n";
  //int pos_ = p_.second;
  if (!invalid_) {
    //std::cout << "IS NOT INVALID\n";
    //if (s_[pos_+Kmer::k] == 0) {
    //                                                                            DID SOMETHING HERE, added -1
    //std::cout <<" NEW ++\n";
    //if (s_[p_.second+total_length_-1] == 0) {
    if (s_[p_.second+min_len_+min_window_-1] == 0) {  
      //std::cout << "IS TOO SMALL\n";
      invalid_ = true;
      return *this;
    } else {
      //std::cout << "IS GOOD, FIND NEXT\n";
      while(!find_next_III()){
        p_.second += 1;
        //std::cout << p_.second << "\n";
      }
      //while(!find_next_II(pos_,pos_+total_length_-1,true)){pos_ = p_.second;}
      //find_next(pos_,pos_+Kmer::k-1,true);
      return *this;
    }
  }
  return *this;
}



// use:  iter++;
// pre:  
// post: iter has been incremented by one
KmerIterator KmerIterator::operator++(int) {
  //std::cout << "ITERATOR ++ V1 OPERATION STARTS\n";
  KmerIterator tmp(*this); 
  operator++(); 
  return tmp;
}


// use:  val = (a == b);
// pre:   
// post: (val == true) if a and b are both exhausted
//       OR a and b are in the same location of the same string.
//       (val == false) otherwise.
bool KmerIterator::operator==(const KmerIterator& o) {
  if (invalid_  || o.invalid_) {
    return invalid_ && o.invalid_;
  } else {
    return (s_ == o.s_) && (p_.second == o.p_.second);
  }
}


// use:  p = *iter;
// pre:   
// post: p is NULL or a pair of Kmer and int
std::pair<Kmer, int>& KmerIterator::operator*() {
  return p_;
}


// use:  example 1: km = iter->first; 
//       example 2:  i = iter->second;
// pre:  *iter is not NULL
// post: km will be (*iter).first, i will be (*iter).second
std::pair<Kmer, int>* KmerIterator::operator->() {
  return &(operator*());
}


// use:  iter.raise(km, rep);
// post: iter has been incremented by one
//       if iter is not invalid, km is iter->first and rep is km.rep()
void KmerIterator::raise(Kmer &km, Kmer &rep) {
  operator++();
  if (!invalid_) {
    km = p_.first;
    rep = km.rep();
  }
}

// use:  find_next(i,j, last_valid); 
// pre:  
// post: *iter is either invalid or is a pair of:
//       1) the next valid kmer in the string that does not have any 'N'
//       2) the location of that kmer in the string
void KmerIterator::find_next(size_t i, size_t j, bool last_valid)
{
  ++i;
  ++j;

  while (s_[j] != 0) 
  {
    char c = s_[j];
    if (c == 'A' || c == 'C' || c == 'G' || c == 'T') 
    {
      if (last_valid) 
      {
        p_.first = p_.first.forwardBase(c);
        break; // default case, 
      } 
      else 
      {
        if (i + Kmer::k - 1 == j) 
        {
          p_.first = Kmer(s_+i);
          last_valid = true;
          break; // create k-mer from scratch
        }
        else
        {
          ++j;
        }
      }
    } 
    else
    {
      ++j;
      i = j;
      last_valid = false;
    }
  }
  if (i+Kmer::k-1 == j && s_[j] != 0)
  {
    p_.second = i;
  } 
  else
  {
    invalid_ = true;
  }
}

//std::tuple<Kmer, bool> KmerIterator::extract_spaced_kmer()
bool KmerIterator::extract_spaced_kmer()
{
  //std::cout << "SUMIMASEN WINDOW SIZE IS " << window_.size() << std::endl;
  //Kmer spaced_kmer = Kmer();

  std::string forward = "";
  std::string backward = "";
  char const * garbo = "A";

  for(uint e = 0; e < total_length_; e+=1)
  {
    if(care_status_->at(e))
    {
      char c = window_.at(e);
      if (c == 'A'){
        forward = forward + 'A';
        backward = 'T' + backward;
      } else if (c == 'C'){
        forward = forward + 'C';
        backward = 'G' + backward;
      } else if (c == 'G'){
        forward = forward + 'G';
        backward = 'C' + backward;
      } else if (c == 'T'){
        forward = forward + 'T';
        backward = 'A' + backward;
      } else {
        p_.first = Kmer(garbo);
        return false;
        //return std::make_tuple(spaced_kmer, false);
      }
    }
  }

  char canonical[forward.length()+1];

  if (forward.compare(backward) < 0){
    //std::cout << forward << std::endl;
    //canonical = &forward[0];
    strcpy(canonical, forward.c_str());
  } else {
    strcpy(canonical, backward.c_str());
    //std::cout << backward << std::endl;
    //canonical = &backward[0];
  }

  p_.first = Kmer(canonical);

  //std::cout << canonical << std::endl;

  //std::cout << "WE ARE FINE" << std::endl;
  //return std::make_tuple(spaced_kmer, true);
  return true;
}



void KmerIterator::minimizer_check(int index)
{
  int nis = minimizers_.size();
  int nif = minimizers_faster_.size();
  int nim = minimizers_faster_super_.size();

  if (nif != nis || nis != nim)
  {
    std::cout << "**ERROR** DIFFERENT NUMBER OF MINIMIZERS CALCULATED IN READ " <<  index <<  " :: " << nis << " :: "<< nif << " :: " << nim << " ::\n";

    return;
  }

  for (int i = 0; i < nis; i+= 1)
  {
    if (std::get<0>(minimizers_[i][0]) != std::get<0>(minimizers_faster_[i][0]) || std::get<0>(minimizers_[i][0]) != std::get<0>(minimizers_faster_super_[i][0]))
    {
      std::cout << "**ERROR** DIFFERENT MINIMIZERS IN WINDOW " << i <<  " :: " << std::get<0>(minimizers_[i][0]) << " :: " << std::get<0>(minimizers_faster_[i][0]) << " :: " << std::get<0>(minimizers_faster_super_[i][0]) << " ::\n";

      return;
    }
    
  }
  //std::cout << "MINIMIZERS FULL OK\n";
  return;
}


/*
void KmerIterator::find_minimizers_faster()
{

  //std::cout << "FASTER MINIMIZERS\n";
  // Get read length
  int read_pos = 0;
  while (s_[read_pos] != 0){read_pos+=1;}

  read_len_ = read_pos;

  if (read_len_ < min_len_ + min_window_ - 1){return;}

  // Window statistics
  //uint64_t current_min = std::numeric_limits<uint64_t>::max();
  uint64_t minimum_in_window = std::numeric_limits<uint64_t>::max(); // keeps track of the minimum k-mer in the current window
  int current_min_position = 0;
  int current_orientation = 0;  
  int current_position = 0;
  bool minimizer_found = false;
  //std::queue<std::tuple<int,int,int,uint64_t>> minimizer_queue; // queue for window k-mers
  
  // Initialize minimizer vector
  minimizers_faster_.clear();

  // First, check the first window
  for (int i = 0; i < min_window_; i+=1)
  {
    current_position = i;

    uint64_t fi = string2int64(i, 0); // forward strand
    uint64_t ri = string2int64(i, 1); // reverse strand

    uint64_t f = minimizer_hash(fi); // forward strand
    uint64_t r = minimizer_hash(ri); // reverse strand

    if (f == r){continue;}// if ambiguous strand we skip

    if (f < r && f <= minimum_in_window) 
    {
      minimum_in_window = f;
      current_min_position = current_position;
      current_orientation = 0;
      minimizer_found = true;
    }
    if (r < f && r <= minimum_in_window) 
    {  
      minimum_in_window = r;
      current_min_position = current_position;
      current_orientation = 1;
      minimizer_found = true;
    }
  }
  if (minimizer_found)
  {
    minimizers_faster_[0].push_back(std::make_tuple(minimum_in_window,current_min_position,current_orientation));
  }
  
  // Then check the rest
  int window_start = 1;

  while (window_start +  min_window_ - 1 + min_len_ - 1 < read_len_)
  {
    current_position = window_start + min_window_ -1;

    uint64_t fi = string2int64(current_position, 0); // forward strand
    uint64_t ri = string2int64(current_position, 1); // reverse strand

    uint64_t f = minimizer_hash(fi); // forward strand
    uint64_t r = minimizer_hash(ri); // reverse strand

    if (f < r && f <= minimum_in_window) 
    {
      minimum_in_window = f;
      current_min_position = current_position;
      current_orientation = 0;
      minimizer_found = true;
    }
    if (r < f && r <= minimum_in_window) 
    {  
      minimum_in_window = r;
      current_min_position = current_position;
      current_orientation = 1;
      minimizer_found = true;
    }

    if (minimizer_found && current_min_position >= window_start)
    {
      minimizers_faster_[window_start].push_back(std::make_tuple(minimum_in_window,current_min_position,current_orientation));
      //std::cout << "NO NEED TO RE-CALC\n";
    }

    if (current_min_position < window_start) // Calculate window minimu again
    {
      //std::cout << "RE-CALC ALL\n";
      minimum_in_window = std::numeric_limits<uint64_t>::max(); // keeps track of the minimum k-mer in the current window
      current_min_position = window_start;
      current_orientation = 0;
      minimizer_found = false;

      // First, check the first window
      for (int j = 0; j < min_window_; j+=1)
      {
        int i = window_start + j;

        current_position = i;

        uint64_t fi = string2int64(i, 0); // forward strand
        uint64_t ri = string2int64(i, 1); // reverse strand

        uint64_t f = minimizer_hash(fi); // forward strand
        uint64_t r = minimizer_hash(ri); // reverse strand

        if (f == r){continue;}// if ambiguous strand we skip

        if (f < r && f <= minimum_in_window) 
        {
          minimum_in_window = f;
          current_min_position = current_position;
          current_orientation = 0;
          minimizer_found = true;
        }
        if (r < f && r <= minimum_in_window) 
        {  
          minimum_in_window = r;
          current_min_position = current_position;
          current_orientation = 1;
          minimizer_found = true;
        }
      }
      if (minimizer_found)
      {
        minimizers_faster_[window_start].push_back(std::make_tuple(minimum_in_window,current_position,current_orientation));
      }
    }

    window_start += 1;
  }
}
*/


// 10% faster than regular fast IF works correctly
void KmerIterator::find_minimizers_faster_super()
{

  // GET READ LENGTH
  int read_pos = 0;
  while (s_[read_pos] != 0){read_pos+=1;}
  read_len_ = read_pos;

  // RETURN IF IS READ TOO SHORT
  if (read_len_ < min_len_ + min_window_ - 1){return;}

  // SOME VARIABLES TO TRACK CURENT MINIMIZERS
  uint64_t minimum_in_window = std::numeric_limits<uint64_t>::max(); // keeps track of the minimum k-mer in the current window
  uint64_t minimum_in_window_reverse = std::numeric_limits<uint64_t>::max();
  int current_min_position = 0;
  int current_min_position_reverse = 0;
  int current_orientation = 0;
  int current_orientation_reverse = 1;  
  int current_position = 0;
  bool minimizer_found = false;
  bool minimizer_found_reverse = false;

  // QUEUES THAT CONTAIN ALL K-MERS IN A WINDOW (AND ITS REVERSE)
  std::deque<std::tuple<int,int,int,uint64_t>> kmer_queue; // queue for window k-mers
  std::deque<std::tuple<int,int,int,uint64_t>> kmer_queue_reverse; // queue for window k-mers
  
  // CLEAR MINIMIZER VECTORS (THESE END UP HAVING THE FINAL MINIMIZERS FOR THE WINDOWS)
  minimizers_faster_super_.clear();
  minimizers_faster_super_reverse_.clear();

  // FIRST, WE FIND THE MINIMIZER FOR THE FIRST WINDOW
  for (int i = 0; i < min_window_; i+=1)
  {
    current_position = i;

    uint64_t f = string2int64(i, 0); // forward strand
    uint64_t r = string2int64(i, 1); // reverse strand

    kmer_queue.push_back(std::make_tuple(1,0,current_position,f));
    kmer_queue_reverse.push_back(std::make_tuple(1,1,current_position,r));

    if (minimizer_hash(f) <= minimizer_hash(minimum_in_window)) 
    {
      minimum_in_window = f;
      current_min_position = current_position;
      current_orientation = 0;
      minimizer_found = true;
    }
    if (minimizer_hash(r) <= minimizer_hash(minimum_in_window_reverse)) 
    {  
      minimum_in_window_reverse = r;
      current_min_position_reverse = current_position;
      current_orientation_reverse = 1;
      minimizer_found_reverse = true;
    }
  }
  if (minimizer_found){minimizers_faster_super_[0].push_back(std::make_tuple(minimum_in_window,current_min_position,current_orientation));}
  if (minimizer_found_reverse){minimizers_faster_super_reverse_[0].push_back(std::make_tuple(minimum_in_window_reverse,current_min_position_reverse,current_orientation_reverse));}
  
  // NEXT, WE FIND THE MINIMIZERS FOR ALL THE FOLLOWING WINDOWS (SHIFTING THE CURRENT WINDOW CHARACTER BY CHARACTER)
  int window_start = 1;

  while (window_start +  min_window_ - 1 + min_len_ - 1 < read_len_)
  {
    current_position = window_start + min_window_ -1;

    uint64_t f = string2int64(current_position, 0); // forward strand
    uint64_t r = string2int64(current_position, 1); // reverse strand

    kmer_queue.push_back(std::make_tuple(1,0,current_position,f));
    kmer_queue_reverse.push_back(std::make_tuple(1,1,current_position,r));
    
    kmer_queue.pop_front();
    kmer_queue_reverse.pop_front();

    if (minimizer_hash(f) <= minimizer_hash(minimum_in_window)) 
    {
      minimum_in_window = f;
      current_min_position = current_position;
      current_orientation = 0;
      minimizer_found = true;
    }
    if (minimizer_hash(r) <= minimizer_hash(minimum_in_window_reverse))
    {  
      minimum_in_window_reverse = r;
      current_min_position_reverse = current_position;
      current_orientation_reverse = 1;
      minimizer_found_reverse = true;
    }

    if (current_min_position < window_start){minimizer_found = false;}
    if (current_min_position_reverse < window_start){minimizer_found_reverse = false;}

    if (minimizer_found){minimizers_faster_super_[window_start].push_back(std::make_tuple(minimum_in_window, current_min_position, current_orientation));}
    if (minimizer_found_reverse){minimizers_faster_super_reverse_[window_start].push_back(std::make_tuple(minimum_in_window_reverse, current_min_position_reverse, current_orientation_reverse));}

    // CALCULATE WINDOW MINIMIZER FOR FORWARD STRAND FROM SCRATCH BECAUSE THE MINIMIZER DROPPED OUT OF THE WINDOW
    if (current_min_position < window_start)
    {
      minimum_in_window = std::numeric_limits<uint64_t>::max(); // keeps track of the minimum k-mer in the current window
      current_min_position = 0;
      current_orientation = 0;
      minimizer_found = false;

      for (int apo = 0; apo < kmer_queue.size(); apo+=1)
      {
        auto entry = kmer_queue.at(apo);
        if (std::get<0>(entry) == 0){continue;}
        if (minimizer_hash(std::get<3>(entry)) <= minimizer_hash(minimum_in_window))
        {
          current_orientation = std::get<1>(entry);
          current_min_position = std::get<2>(entry);
          minimum_in_window = std::get<3>(entry);
          minimizer_found = true;
        }
      }
      if (minimizer_found){minimizers_faster_super_[window_start].push_back(std::make_tuple(minimum_in_window, current_min_position, current_orientation));}
    }

    // CALCULATE WINDOW MINIMIZER FOR REVERSE STRAND FROM SCRATCH BECAUSE THE MINIMIZER DROPPED OUT OF THE WINDOW
    if (current_min_position_reverse < window_start)
    {
      minimum_in_window_reverse = std::numeric_limits<uint64_t>::max(); // keeps track of the minimum k-mer in the current window
      current_min_position_reverse = 0;
      current_orientation_reverse = 0;
      minimizer_found_reverse = false;

      for (int apo = 0; apo < kmer_queue_reverse.size(); apo+=1)
      {
        auto entry = kmer_queue_reverse.at(apo);
        if (std::get<0>(entry) == 0){continue;}
        if (minimizer_hash(std::get<3>(entry)) <= minimizer_hash(minimum_in_window_reverse))
        {
          current_orientation_reverse = std::get<1>(entry);
          current_min_position_reverse = std::get<2>(entry);
          minimum_in_window_reverse = std::get<3>(entry);
          minimizer_found_reverse = true;
        }
      }
     if (minimizer_found_reverse){minimizers_faster_super_reverse_[window_start].push_back(std::make_tuple(minimum_in_window_reverse, current_min_position_reverse, current_orientation_reverse));}
    }

    // MOVE TO NEXT WINDOW
    window_start += 1;
  }
  // MINIMIZER ARE NOW CALCULATED
}


/*
void KmerIterator::find_minimizers()
{
  // Get read length
  imfine = true;
  int read_pos = 0;
  while (s_[read_pos] != 0){read_pos+=1;}

  read_len_ = read_pos;

  // Initializr queue for minimizers (not in use yet)
  //std::queue<std::tuple<int,int,int>> minimizer_queue;

  // Initialize minimizer vector
  minimizers_.clear();

  for (int i = 0; i < read_len_ - min_len_ - min_window_ + 2; i+=1) // for every minimizer window position
  {
    //std::cout << "MINIMIZER SEARCH POSITION: " << i << "\n";

    uint64_t m = std::numeric_limits<uint64_t>::max(); // max hash value

    for (int j = 0; j < min_window_; j+= 1) // for every k-mer in the minimizer window
    {
      uint64_t fi = string2int64(i+j, 0); // forward strand
      uint64_t ri = string2int64(i+j, 1); // reverse strand

      uint64_t f = minimizer_hash(fi); // forward strand
      uint64_t r = minimizer_hash(ri); // reverse strand

      //std::cout << "F VALUE " << fi << " F HASH " << f << "\n";
      //std::cout << "R VALUE " << ri << " R HASH " << r << "\n";

      if (f != r) // if ambiguous strand we skip
      {
        m = std::min(m, std::min(f,r));
      }
    }

    
    for (int j = 0; j < min_window_; j+= 1) // for every k-mer in the minimizer window
    {
      uint64_t fi = string2int64(i+j, 0); // forward strand
      uint64_t ri = string2int64(i+j, 1); // reverse strand

      uint64_t f = minimizer_hash(fi); // forward strand
      uint64_t r = minimizer_hash(ri); // reverse strand

      if ((f == m) & (f < r))
      {
        minimizers_[i].push_back(std::make_tuple(f,i+j,0));
      }
      else if ((r == m) & (r < f))
      {
        minimizers_[i].push_back(std::make_tuple(r,i+j,1));
      }
    }
  }

  /*
  std::cout << "MINIMIZERS CALCULATED, THE POSITIONS ARE\n";

  for(auto const& imap: minimizers_)
  {
    std::cout << imap.first << "\n";
    //vints.push_back(imap.first);
  }
  
  //exit(10);
}
*/



/*
    Finds next k-mer using the window minimizers
*/

bool KmerIterator::find_next_III()
{
  std::vector<std::tuple<int, int> > fixed_streak_starts; //orientation (0 forward, 1 reverse) AND position
  std::vector<std::tuple<int, int> > fixed_streak_starts_reverse;

  int position_tmp = -1;
  int orientation_tmp = -1;
  int position_tmp_reverse = -1;
  int orientation_tmp_reverse = -1;

  for (int i = 0; i < min_num_; i+=1) // The number of minimizers needed for a spaced k-mer
  { 
    //int min_pos = p_.second + i*(min_len_+min_window_-1); // position of the next minimizer
    int min_pos = p_.second + i*(min_window_offset_); // position of the next minimizer

    if (min_pos >= read_len_)
    {
      invalid_ = true;
      return true;
    }

    if (minimizers_faster_super_.count(min_pos) < 1 && minimizers_faster_super_reverse_.count(min_pos) < 1){return false;}
    
    if (minimizers_faster_super_.count(min_pos) > 0)
    {
      position_tmp = std::get<1>(minimizers_faster_super_[min_pos][0]);
      orientation_tmp = std::get<2>(minimizers_faster_super_[min_pos][0]);
      fixed_streak_starts.push_back(std::make_tuple(orientation_tmp, position_tmp));
    }

    if (minimizers_faster_super_reverse_.count(min_pos) > 0)
    {
      position_tmp_reverse = std::get<1>(minimizers_faster_super_reverse_[min_pos][0]);
      orientation_tmp_reverse = std::get<2>(minimizers_faster_super_reverse_[min_pos][0]);
      fixed_streak_starts_reverse.push_back(std::make_tuple(orientation_tmp_reverse, position_tmp_reverse));
    } 
  }


  bool forward_is_valid = false;
  if (fixed_streak_starts.size() == min_num_){forward_is_valid = true;}

  bool reverse_is_valid = false;
  if (fixed_streak_starts_reverse.size() == min_num_){reverse_is_valid = true;}



  // Build FORWARD spaced k-mer by concatenating fixed streaks
  char minstring[min_len_*min_num_+1];  
  int fmwp = 0; // forward minimizer window position
  int mcp = 0; // minimizer character position

  for (auto mpot : fixed_streak_starts)
  {
    int mo = std::get<0>(mpot);
    int mp = std::get<1>(mpot);
    //std::cout << "Orientation " << mo << " Position " << mp << "\n";
    // FORWARD
    for (int x = 0; x < min_len_; x+=1)
    {
      minstring[fmwp+mcp] = s_[mp+x];
      mcp+=1;
    }
    mcp=0;
    fmwp = fmwp + min_len_;
  }
  minstring[min_len_*min_num_] = '\0';
  //std::cout << minstring << "\n";
  Kmer next_kmer_candidate1 = Kmer(minstring);
  std::string ks = next_kmer_candidate1.toString();


  // Build REVERSE spaced k-mer by concatenating fixed streaks
  char minstring_reverse[min_len_*min_num_+1];  
  int fmwp_reverse = 0; // forward minimizer window position
  int mcp_reverse = 0; // minimizer character position

  for (auto mpot : fixed_streak_starts_reverse)
  {
    int mo = std::get<0>(mpot);
    int mp = std::get<1>(mpot);
    //std::cout << "Orientation " << mo << " Position " << mp << "\n";
    // FORWARD
    for (int x = 0; x < min_len_; x+=1)
    {
      minstring_reverse[fmwp_reverse+mcp_reverse] = s_[mp+x];
      mcp_reverse+=1;
    }
    mcp_reverse=0;
    fmwp_reverse = fmwp_reverse + min_len_;
  }
  minstring_reverse[min_len_*min_num_] = '\0';
  //std::cout << minstring_reverse << "\n";
  Kmer next_kmer_candidate1_reverse_flipped = Kmer(minstring_reverse);                                                                    // SUS!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!11
  Kmer next_kmer_candidate1_reverse = next_kmer_candidate1_reverse_flipped.twin_MOD();
  std::string ks_reverse = next_kmer_candidate1_reverse.toString();


  Kmer true_minimizer_kmer;

  if (forward_is_valid && reverse_is_valid)
  {
    if (next_kmer_candidate1.toString() <= next_kmer_candidate1_reverse.toString()){
      true_minimizer_kmer = next_kmer_candidate1;
    } else {
      true_minimizer_kmer = next_kmer_candidate1_reverse;
    }

  }
  else if (forward_is_valid && !reverse_is_valid)
  {
    true_minimizer_kmer = next_kmer_candidate1;
  }
  else if (!forward_is_valid && reverse_is_valid)
  {
    true_minimizer_kmer = next_kmer_candidate1_reverse;
  }
  else
  {
    return false;
  }

  std::string tks = true_minimizer_kmer.toString();

  // BLOCK BELOW: Discard k-mer if it was found in the PREVIOUS position of the read

  if (previous_minimizer_kmer_ == tks)
  {
    return false;
  }
  else if (found_minimizer_kmers_.count(tks) != 0)
  {
    return false;
  }
  else
  {
    found_minimizer_kmers_[tks] = 1;
    previous_minimizer_kmer_ = tks;
    p_.first = true_minimizer_kmer;
    return true;
  }
}


char KmerIterator::reverse_complement_character(char c)
{
  if (c == 'A'){return 'T';}
  else if (c == 'T'){return 'A';}
  else if (c == 'C'){return 'G';}
  else if (c == 'G'){return 'C';}
  else {std::cout << "** ERROR ** Invalid character, does not have reverse complement\n"; exit(999);}
  return 'Z';
}


/*
    Turns a string of characters to 64 bit integer representation (can handle only strings of length 32 or less)
*/

uint64_t KmerIterator::string2int64(int start, int strand)
{
  uint8_t a = 0;
  uint8_t c = 1;
  uint8_t g = 2;
  uint8_t t = 3;

  uint64_t str_as_int = 0;

  if (strand == 0) // forward strand
  {
    for (int i = start; i < start+min_len_; i+=1)
    {
      str_as_int = str_as_int << 2;
      if (s_[i] == 'A'){str_as_int = str_as_int | a;}
      else if (s_[i] == 'C'){str_as_int = str_as_int | c;}
      else if (s_[i] == 'G'){str_as_int = str_as_int | g;}
      else if (s_[i] == 'T'){str_as_int = str_as_int | t;}
      else {
        imfine=false;
        return std::numeric_limits<uint64_t>::max();
        //std::cout << "ERROR CHARACTER FORWARD " << s_[i] << " \n";
      }
      //std::cout << "KMER NOW: " << str_as_int << "\n";
    }
    return str_as_int;
  }
  else if (strand == 1) // reverse strand
  {
    int back = 0;
    for (int j = 0; j < min_len_; j+=1)
    {
      int i = start + min_len_ - 1 - j;
      str_as_int = str_as_int << 2;
      if (s_[i] == 'A'){str_as_int = str_as_int | t;}
      else if (s_[i] == 'C'){str_as_int = str_as_int | g;}
      else if (s_[i] == 'G'){str_as_int = str_as_int | c;}
      else if (s_[i] == 'T'){str_as_int = str_as_int | a;}
      else {
        imfine=false;
        return std::numeric_limits<uint64_t>::max();
        //std::cout << "ERROR CHARACTER REVERSE " << s_[i] << " \n";
      }
      //std::cout << "KMER NOW: " << str_as_int << "\n";
    }
    return str_as_int;
  }
  std::cout << "RETURNING MAX, TRIED TO USE STRAND THAT IS NOT FORWARD OR BACKWARD\n";
  return std::numeric_limits<uint64_t>::max();
}


/*
    COPIED FROM https://naml.us/post/inverse-of-a-hash-function/

    ===== REMEMBER TO CITE =====

*/
uint64_t KmerIterator::minimizer_hash(uint64_t key) {
  

  if (key ==  std::numeric_limits<uint64_t>::max())
  {
    //std::cout << "### Minimizer hash ERROR ###\n";                                                                        // DISABLED NOW BECAUSE WE GO HERE
    return std::numeric_limits<uint64_t>::max();
  }

  // FOR TESTS
  //return key;

  key = (~key) + (key << 21); // key = (key << 21) - key - 1;
  key = key ^ (key >> 24);
  key = (key + (key << 3)) + (key << 8); // key * 265
  key = key ^ (key >> 14);
  key = (key + (key << 2)) + (key << 4); // key * 21
  key = key ^ (key >> 28);
  key = key + (key << 31);
  if (key == 0){std::cout << "ZERO KEY\n";}
  return key;
}




// DOES NOT WORK WHEN INITIALIZATION FAILS ON THE FIRST POSITION

// use:  find_next(i,j, last_valid); 
// pre:  
// post: *iter is either invalid or is a pair of:
//       1) the next valid kmer in the string that does not have any 'N'
//       2) the location of that kmer in the string
//
//
// NEW MODIFIED VERSION
//
bool KmerIterator::find_next_II(size_t i, size_t j, bool last_valid)
{
  ++i; // Start position
  ++j; // End position

  

  // If string does not have enough characters left
  //if (sizeof(s_)/sizeof(char) - i < total_length_)
  //{
  //  invalid_ = true;
  //  return;
  //}

  //std::cout << "HALOO" << std::endl;

  // If we do not have a long k-mer ready, build it
  if(!last_valid)
  {
    //std::cout << "ORA ORA" << std::endl;
    for(uint si = i; si < i+total_length_; si+=1)
    {
      if (s_[si] == 0){invalid_=true;return true;}
      //std::cout << "HALOOOOOO" << std::endl;
      char c = s_[si];
      if (c == 'A' || c == 'C' || c == 'G' || c == 'T')
      {
        window_.push_back(c);
      }
      else
      {
        window_.push_back('N');
      }
    }
  }
  // Otherwise add just one character

  else
  {
    //std::cout << "ARA ARA" << std::endl;
    if (s_[i+total_length_-1] == 0){invalid_=true; return true;}
    //std::cout << "HALOOOOOO" << std::endl;
    char c = s_[i+total_length_-1];
    if (c == 'A' || c == 'C' || c == 'G' || c == 'T')
    {
      window_.push_back(c);
    }
    else
    {
      window_.push_back('N');
    }
    window_.pop_front();
  }

  bool kmer_fine;

  kmer_fine = extract_spaced_kmer(); 

  if(kmer_fine)
  {
    //p_.first = spaced_kmer;
    p_.second = i;
    return true;
  }
  else
  {
    //std::cout << "SONNA\n";
    p_.second = i;
    return false;
  }
}
