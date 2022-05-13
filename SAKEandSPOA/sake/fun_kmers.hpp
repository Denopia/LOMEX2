#ifndef FUN_KMERS_H
#define FUN_KMERS_H

#include <iostream>
#include <fstream>
#include <vector>
using namespace std;





/*
	Map nucleotides to integers:
		C = 0 (00 as two bits)
		A = 1 (01 as two bits)
		T = 2 (10 as two bits)
		G = 3 (11 as two bits)
*/

//enum NUC_MAP {C, A, T, G};

enum NUC_MAP {A, C, G, T};





/*
	Given an integer, return the corresponding
	nucleotide character
*/
char map_int2nuc(uint8_t nuc);


/*
	Given a character, return the corresponding
	nucleotide integer
*/
uint8_t map_nuc2int(char nuc);


/*
	Given a 128 bit integer, return the
	corresponding nucleotide string
*/
string map_int2str(__uint128_t seq, uint64_t len);


/*
	Given a string, return the corresponding
	nucleotide sequence as 128 bit integer
*/
__uint128_t map_str2int(std::string seq);




/*
	Given a string, return the corresponding
	nucleotide sequence as 128 bit integer
*/
void map_str2intR(std::string & seq, __uint128_t & bin_seq);


/*
	Given a nucleotide as an integer,
	return its reverse complement as an integer
*/
int reverse_complement_nucint(int nuc);


/*
	Given a nuclotide as a character,
	return its reverse complement as a character
*/

char reverse_complement_nucchar(char nuc);


/*
	Given a nucleotide sequence as a 128 bit integer
	and the number of the nucleotides, return its 
	reverse complement as 128 bit integer 
*/
__uint128_t reverse_complement_seqint(__uint128_t seq, uint64_t len);


/*
	Given a nuclotide sequence as a string,
	return its reverse complement as a string
*/
std::string reverse_complement_seqstr(string seq);


/*
	Given two sequences, determines if the first one is
	bigger according to the mapped nuclotide values.
	(Used to compare k-mer and its reverse complement)
*/
bool compare_seqs(__uint128_t seq1, __uint128_t seq2);


/*
	Given two sequences, determines if the first one is
	bigger according to the mapped nuclotide values.
	(Used to compare k-mer and its reverse complement)
*/
bool compare_seqs_string(std::string seq1, std::string seq2);


/*
	Given a 64 bit integer x and another integer y,
	checks if the y rightmost bits in x are zeros.
*/
bool is_zero(uint64_t x, int y);


/*
	Determines if a set of sequences contains characters
	other than C, A, T or G.

	Input: 
		Array of 64 bit integers, one for every block 
		of continuous fixed character sequences. One bit
		in the integers corresponds to a single nucleotide
		in the corresponding block starting from the
		rightmost bit. If the bit is 1, that character is
		not in the standard nucleotides C,A,T,G.

		Vector of integers, which tells how many characters
		each fixed character block has.

		Integer, which tells how many fixed character blocks
		there are in total.

	Output:
		True or false. If any of the fixed blocks contains
		non-standard nucleotides, return false. Otherwise
		all characters are C, A, T or G, and true is returned.
*/
bool no_Ns_present(uint64_t nbs[], vector<int> &lengths, int blocks);


/*
	Given a spaced seed pattern as a string, return 
	the corresponding boolean vector. Returns a tuple,
	which contains the boolean vector, total length of
	the pattern, number of fixed positions, the number
	of continuous fixed and don't care blocks, and
	the starting positions of these blocks.
*/
tuple<vector<bool>, int, int, vector<int>, vector<int>, vector<int> > interpret_spaced_seed_pattern(string pattern);



/*
	Extracts spaced k-mer from a long k-mer using the given spaced k-mer pattern
*/
std::string extract_spaced_kmer(std::string long_kmer, std::vector<bool> & is_fixed_character);




/*
	Given a string, return the corresponding
	nucleotide sequence as 128 bit integer
*/
uint64_t map_str2int_small(std::string seq);



std::string map_int2str_small(uint64_t seq, int len);

uint64_t reverse_complement_seqint_small(uint64_t seq, int len);

/*
	Given an integer, return its hash
*/
uint64_t minimizer_hash(uint64_t key);


/*
	Minimizer searcher
*/
std::tuple<uint64_t, int, int, bool> find_minimizer_from_scratch(std::vector<uint64_t> & current_read_kmers_as_ints, int window_start_pos, int window_size, std::vector<bool> & has_n);


/* Give random nucleotide character */
char random_nuc();

#endif
