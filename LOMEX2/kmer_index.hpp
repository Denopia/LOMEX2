#ifndef KMER_INDEX_H
#define KMER_INDEX_H


#include <iostream>
#include <fstream>
#include <tuple>
#include <vector>
#include <deque>
#include <map>
#include <unordered_set>



class KMerIndex
{

private:
	int kmers_n;
	int kmer_len;
	int occurrences_n;
	std::string kmers_path;
	std::map<uint64_t, std::vector<int> > kmers2occurrences;
	std::unordered_set<uint64_t> relevant_kmers;


public:

	void initialize_variables(std::string input_kmers_path, int input_kmer_len);

	void initialize_index(int minimizer_set_occ_threshold);

	void add_occurrence(uint64_t kmer, int read_id, int position, int length, int strand);

	uint64_t canonical_kmer(uint64_t kmer);

	void print_index();

	bool is_this_relevant_kmer(uint64_t kmer);

	int get_kmers_n(){return kmers_n;}

	std::unordered_set<uint64_t>& get_kmer_set(){return relevant_kmers;}

	std::vector<int>& get_kmer_info(uint64_t kmer){return kmers2occurrences[kmer];} 

};

std::tuple<std::string, int> interpret_minimizer_set_line(std::string one_line);


#endif