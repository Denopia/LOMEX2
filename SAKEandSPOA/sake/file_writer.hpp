#ifndef FILE_WRITER_H
#define FILE_WRITER_H

#include <iostream>
#include <fstream>
#include <tuple>
#include <vector>
#include <deque>
#include <map>
#include <random>

using namespace std;


/*
	Write locations into a file
*/
std::string write_locations(map<__uint128_t, vector<int> > & kmer2positions, std::string work_dir, int file_number, int fixed_length);


/*
	Write occurrences into a file

	SUPPORT ONLY FOR k-MERS UP TO k = 255*4 = 1020 

*/
std::string write_occurrences_binary(map<__uint128_t, vector<std::string> > & kmer2occurrences, std::string work_dir, int file_number, int fixed_length, int total_length);


/*
	Write occurrences into a file

	SUPPORT ONLY FOR k-MERS UP TO k = 255*4 = 1020 

*/
std::string write_occurrences_binary_dupesnt(map<__uint128_t, vector<std::string> > & kmer2occurrences, std::string work_dir, int file_number, int fixed_length, int total_length);


/*
	Write regular k-mers in a file as binary

*/
tuple<std::string, uint64_t> write_regular_kmers_ready_binary(map<__uint128_t, vector<uint8_t> > & spaced2regular, std::string work_dir, int file_number, int fixed_length, int total_length, int thread);


/*
	Write regular k-mers in a file as binary (for debugging)

*/
tuple<std::string, uint64_t> write_regular_kmers_ready_binary_dupesnt_DEBUG(map<__uint128_t, vector<uint8_t> > & spaced2regular, std::string work_dir, int file_number, int fixed_length, int total_length, int thread, vector<bool> & character_status);


/*
	Write regular k-mers in a file as binary

*/
tuple<std::string, uint64_t> write_regular_kmers_ready_binary_dupesnt(map<__uint128_t, vector<uint8_t> > & spaced2regular, std::string work_dir, int file_number, int fixed_length, int total_length, int thread, vector<bool> & character_status);



/*
	Put regular k-mer into a buffer as bytes forward 

*/
void put_kmer_in_buffer_forward(std::vector<char> & read_vector, int ri, int total_length, map<__uint128_t, vector<uint8_t> > & spaced2regular, __uint128_t stored_read_spaced_kmer);


/*
	Put regular k-mer into a buffer as bytes backward

*/
void put_kmer_in_buffer_backward(std::vector<char> & read_vector, int ri, int total_length, map<__uint128_t, vector<uint8_t> > & spaced2regular, __uint128_t stored_read_spaced_kmer);


/*
	Put regular k-mer into a buffer (forward)

*/
void put_kmer_in_buffer(std::deque<uint8_t> & push_kmer, map<__uint128_t, vector<uint8_t> > & spaced2regular, __uint128_t stored_read_spaced_kmer);

/*
	Format the file number
*/
std::string format_number(int number);


/*
	Maps 4 nucleotide characters into a byte
*/
uint8_t map_4nucs2byte(char nucs[4]);


/*
	Gives a random nucleotide character
*/
char random_nucleotide_character();


/*
	Gives a random nucleotide integer
*/
uint8_t random_nucleotide_int();


/*
	Map a byte into four nucleotides
*/
std::string map_byte2fournucs(uint8_t byte);



// New experimental stuff
void put_kmer_in_buffer_REE(std::deque<uint8_t> & long_kmer, map<__uint128_t, vector<uint8_t> > & spaced2regular, __uint128_t stored_read_spaced_kmer);

// New experimental stuff
void update_last4(std::deque< std::deque<uint8_t> > & last4longkmers, std::deque<uint8_t> & long_kmer_as_bits, int last_amigos);


#endif