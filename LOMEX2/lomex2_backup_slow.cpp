#include <boost/program_options.hpp>
#include <boost/log/trivial.hpp>
#include <boost/algorithm/string.hpp>
#include <boost/algorithm/string/predicate.hpp>
#include <boost/filesystem.hpp>
#include <boost/range/iterator_range.hpp>
#include <map>
#include <math.h>
#include <iostream>
#include <fstream>
#include <sstream>
#include <bitset>
#include <regex>
#include <string>
#include <tuple>
#include <vector>
#include <deque>
#include <random>
#include <map>
#include <chrono>
#include <limits>
#include <math.h>
#include <cstdlib>
#include "fun_kmers.hpp"
#include "file_writer.hpp"
#include "file_reader.hpp"
#include "kmer_index.hpp"
#include "fun_poa.hpp"


/*
extern "C" 
{
 	#include "lpo.h"
	#include "msa_format.h"
	#include "align_score.h"
}
*/

//namespace chrono = std::chrono; 
namespace po = boost::program_options;
namespace bfs = boost::filesystem;
//using namespace std;
//using namespace boost::filesystem;


int main(int argc, char *argv[])
{

	std::vector<std::string> test_sequences;
	std::vector<std::string> consensus_sequences;
	float bundling_threshold = 0.8;

	test_sequences.push_back("AATAATTTAATAA");
	test_sequences.push_back("AATATTTAATAA");
	test_sequences.push_back("AATAATTTAAATAA");
	test_sequences.push_back("AATAATTAATAA");

	LOMEXPOAGlue lpoag;
	lpoag.initialize_me("blosum80.mat");

	//consensus_sequences = lpoag.poa_consensus_sequences(test_sequences, bundling_threshold);

	//return 0;

	/*
		Initialize variables for user given arguments.
	*/
	std::string reads_path, kmers_path, output_path;
	int window_number, window_size, minimizer_len;
	float placeholder_float;

	/*
		Parse arguments given by the user.
	*/
	po::options_description desc("LOMEXII options");
	desc.add_options()
		("help,h", "Give help")
		("kmers,k", po::value<std::string>(& kmers_path)->default_value("na0"), "Path to the k-mers file")
		("reads,r", po::value<std::string>(& reads_path)->default_value("na1"), "Path to the reads file")
		("output,o", po::value<std::string>(& output_path)->default_value("na2"), "Path to the output file")
		("windows,n", po::value<int>(& window_number)->default_value(1), "Number of windows")
		("window-size,s", po::value<int>(& window_size)->default_value(1), "Size of window")
		("minimizer-len,l", po::value<int>(& minimizer_len)->default_value(1), "Minimizer length")
		("pl-float,f", po::value<float>(& placeholder_float)->default_value(1.1), "Placeholder float");
		
	po::variables_map vm;
	po::store(po::parse_command_line(argc, argv, desc), vm);
	po::notify(vm);


	/*
		If user asks for help, print help message and exit.
	*/
	if (vm.count("help")) {
 		std::cout << "Print help message here" << std::endl;
    	return 0;
	}

	/*
		Later, add a check to make sure arguments are valid.
	*/


	/*
	std::cout << kmers_path << "\n";
	std::cout << reads_path << "\n";
	std::cout << output_path << "\n";
	std::cout << window_number << "\n";
	std::cout << window_size << "\n";
	std::cout << minimizer_len << "\n";
	*/

	//return 0;

	/*
		Start timing the program execution.
	*/

	/*
	int execution_seconds;
	execution_seconds = 0;
	std::chrono::high_resolution_clock::time_point start_time, end_time;
	std::chrono::duration<double> elapsed_time;
	*/


	/*
		Use the code block below to measure execution time.
	*/

	/*
	start_time = std::chrono::high_resolution_clock::now();
	end_time = std::chrono::high_resolution_clock::now();
	elapsed_time = end_time - start_time;
	execution_seconds += (elapsed_time.count());
	*/


	/*
		**** ACTUAL PROGRAM STARTS HERE ****
	*/


	/*
		Map where minimizer-k-mers are stored associated with an integer vector.
		Vector has antries for each occurrence in the reads.
		A single occurrence is a triplet of consecutive integers (r, p, l),
		where r = read id, p = position in the read, l = length of the occurrence
	*/
	//std::map<uint64_t, std::vector<int> > kmer_map;


	int minimizer_kmer_len = window_number * minimizer_len;

	KMerIndex kmer_index;
	kmer_index.initialize_variables(kmers_path, minimizer_kmer_len);
	kmer_index.initialize_index();


	/*
		This is a fastq file reader used to go through the reads.
		(The two last parameters do nothing and can be ignored for now.)
	*/

	std::cout << "Initialize fastq reader\n";

	FastqFileReader reads_file_reader;
	reads_file_reader.initialize_me(reads_path, 0, 0);


	/*

		Now check all reads for minimimizer-k-mers in the kmer index

	*/


	//std::deque<char> current_read; NOT IN USE

	// === Initialize some needed variables ===
	char current_character; // current character is here
	
	std::vector<uint64_t> current_read_kmers_as_ints; //all minimizer length k-mers as integers
	std::vector<uint64_t> current_read_kmer_hahses; // all hashed minimizer length k-mers
	std::vector<uint64_t> current_read_minimizers; // all minimizers
	std::vector<uint64_t> current_read_minimizer_kmers; // all minimizers
	
	uint64_t current_kmer_as_int; // Current k-mer as integer
	int characters_read_from_read; // Simple character counter to keep track of how many characters are read from a read
	uint64_t minimizer_mask = pow(2, 2*minimizer_len)-1; // Mask with ones where minimizer would have characters

	std::string current_kmer_as_string; // k-mer string for debugging


	//----------------------------------------------------------------------
	// Debug stuff
	//
	// 

	std::vector<std::string> all_reads;
	std::string one_read;

	//
	//
	//-------------------------------------------------------------------------

	// === Loop through every read in the input file (MUST BE FASTQ) ===
	while(reads_file_reader.get_reads_left())
	{
		one_read = "";//-------------------------------------------------------------------------------------------------------------------------------------
		// === Read the next read ===
		reads_file_reader.roll_to_next_read();
		// === Clear variables related to the previous read ===
		current_read_kmers_as_ints.clear();
		current_read_kmer_hahses.clear();
		current_read_minimizers.clear();
		current_read_minimizer_kmers.clear();
		
		current_kmer_as_int = 0;
		characters_read_from_read = 0;
		current_kmer_as_string = "";

		std::cout << "This is read number " << reads_file_reader.get_current_read_number() << "/" << reads_file_reader.get_total_reads() << " with length " << reads_file_reader.get_read_length() << "\n";

		// === Loop through all characters of the current read TO FIND ALL K-MERS IN IT===
		while (!reads_file_reader.is_empty())
		{
			// === Read next character and update character counter ===
			current_character = reads_file_reader.pop_next_char();
			one_read = one_read + current_character; //------------------------------------------------------------------------------------------------------
			characters_read_from_read += 1;
			// === Update current k-mer ===
			current_kmer_as_int = current_kmer_as_int << 2;
			current_kmer_as_int = current_kmer_as_int | map_nuc2int(current_character);
			current_kmer_as_int = current_kmer_as_int & minimizer_mask;

			// === Keep track of string for debugging ===
			//current_kmer_as_string = current_kmer_as_string + current_character;
			//while (current_kmer_as_string.length() > minimizer_len){current_kmer_as_string = current_kmer_as_string.substr(1);}

			// === If enough characters are read push the k-mer into the vector ===
			if (characters_read_from_read >= minimizer_len)
			{
				current_read_kmers_as_ints.push_back(current_kmer_as_int);
				//std::cout << "k-mer: " + current_kmer_as_string << " is " << current_kmer_as_int << "\n";
			}
		} // All k-mers in the current read have been found and stored as ints

		all_reads.push_back(one_read);

		//std::cout << "k-mers found as ints successfully\n";
		//std::cout << "There are " << current_read_kmers_as_ints.size() << " k-mers in the read\n\n";

		// === Here we can hash all k-mers ===
		for (int i = 0; i < current_read_kmers_as_ints.size(); i++)
		{
			current_read_kmer_hahses.push_back(current_read_kmers_as_ints[i]);
		}

		//std::cout << "k-mers hashed successfully\n";
		//std::cout << "There are " << current_read_kmer_hahses.size() << " hashed k-mers in the read\n\n";

		// === After all k-mers are in the vector it is time to find all the minimizers ===
		// ================================================================================
		// ================================================================================

		// This is something that can be made faster if we do not check every window separately

		// ================================================================================
		// ================================================================================

		// === Set the minimizer in this variable ===
		uint64_t i_am_the_minimizer;
		int minimizer_position;
		// === Now go through all windows to find the minimizers ===
		for (int i = 0; i+window_size-1 < current_read_kmer_hahses.size(); i++)
		{
			// === Set the minimizer to the maximum value so it will be updated (or at least tied) ===
			i_am_the_minimizer = std::numeric_limits<uint64_t>::max();
			// === Set the minimizer position at the start of the window (it will be updated) ===
			minimizer_position = i;
			// === Now check every k-mer in the window ===
			for (int j = i; j < i+window_size; j++)
			{
				// === If a new minimizer is found, update variables ===
				if (current_read_kmer_hahses[j] <= i_am_the_minimizer)
				{
					i_am_the_minimizer = current_read_kmer_hahses[j];
					minimizer_position = j;
				}
			}
			// === After all k-mers are checked, push the found minimizer in the minimizer vector
			current_read_minimizers.push_back(i_am_the_minimizer);
		}

		//std::cout << "Minimizers found successfully\n";
		//std::cout << "There are " << current_read_minimizers.size() << " minimizers in the vector\n\n";

		// === Now that all minimizers are found, we can find all minimizer-k-mers ===

		// === First clear the minimizer-k-mer vector of previous read ===
		//current_read_minimizer_kmers.clear();

		// Put minimizer-k-mer here
		uint64_t minimizer_kmer;
		// This is the distance between two minimizer windows
		int minimizer_window_distance = minimizer_len + window_size - 1;

		for (int i = 0; i+minimizer_window_distance*(window_number-1) < current_read_minimizers.size(); i++)
		{
			//std::cout << "Minimizer at position: " << i << "\n";
			minimizer_kmer = 0;
			for (int j = 0; j < window_number; j++)
			{
				minimizer_kmer = minimizer_kmer << 2*minimizer_len;
				minimizer_kmer = minimizer_kmer | current_read_minimizers[i+j*minimizer_window_distance];
			}
			current_read_minimizer_kmers.push_back(minimizer_kmer);
		}

		//std::cout << "minimizer-k-mers found successfully\n";
		//std::cout << "There are " << current_read_minimizer_kmers.size() << " minimizer-k-mers in the vector\n\n";


		// === After this we can finally put the minimizer-k-mer occurrences into the minimizer-k-mer index ===
		uint64_t previous_insert_minimizer_kmer = std::numeric_limits<uint64_t>::max();
		uint64_t current_insert_minimizer_kmer;
		int iloop = 0;
		int jloop = 0;
		int same_minimizer_kmer_streak;
		int minimizer_kmer_full_window_size = window_number*(minimizer_len + window_size - 1);
		while(iloop < current_read_minimizer_kmers.size())
		{
			current_insert_minimizer_kmer = current_read_minimizer_kmers.at(iloop);
			same_minimizer_kmer_streak = 0;
			jloop = iloop+1;
			while (jloop < current_read_minimizer_kmers.size())
			{
				if (current_insert_minimizer_kmer != current_read_minimizer_kmers.at(jloop)){break;}
				//std::cout << "Extended minimizer-k-mer\n";
				jloop+=1;
				same_minimizer_kmer_streak+=1;
			}

			if (kmer_index.is_this_relevant_kmer(current_insert_minimizer_kmer))
			{
				//std::cout << "Inserted minimizer-k-mer\n";
				kmer_index.add_occurrence(current_insert_minimizer_kmer, reads_file_reader.get_current_read_number(), iloop, minimizer_kmer_full_window_size + same_minimizer_kmer_streak);
			}
			iloop = jloop;
		}

		std::cout << "\n";
	} // Every read has been read after this




	// Start building consensus k-mers
	for (auto minikmer : kmer_index.get_kmer_set())
	{
		std::cout << "=========================================================\n";
		std::cout << "K-MER " << " OCCURRENCES\n";
		int idk = 0;
		while (idk < kmer_index.get_kmer_info(minikmer).size())
		{
			std::cout << "READ ID: " << kmer_index.get_kmer_info(minikmer)[idk] << "\n";
			idk+=1;
			std::cout << "POSITION IN READ: " << kmer_index.get_kmer_info(minikmer)[idk] << "\n";
			idk+=1;
			std::cout << "LENGTH: " << kmer_index.get_kmer_info(minikmer)[idk] << "\n";
			idk+=1;
			std::cout << "-----------------------------------------------------\n";			
		}

		std::cout << "=========================================================\n";
		std::cout << minikmer << "\n";


	}



	return 0;



	/*
		Before starting to go thru reads, initialize variables to track
		where we are going
	*/

	/*

	int current_read = , current_position;
	int current_position, next_position;
	current_read = -1;

	std::string full_window_string = "";
	int full_window_length = window_number * (minimizer_len + window_size - 1);

	int windows_checked = 0;

	std::string previous_window_minimizer_kmer = "Z";

	while(true)
	{
		//std::cout << "PLS\n";
		std::tie(next_read, next_position, next_character) = fasta_reader.give_next_read_position_and_character();
		

		if (current_read == next_read && current_read != -1){
			full_window_string.push_back(next_character);
			//std::cout << "NANI!\n";
		}
		else {
			full_window_string.clear();
			full_window_string.push_back(next_character);
			//std::cout << "SONNA!\n";
		}

		if (fasta_reader.get_im_finished())
		{
			break;
		}

		current_read = next_read;
		current_position = next_position;
		current_charater = next_character;


		if (full_window_string.length() < full_window_length){continue;}
		while(full_window_string.length() != full_window_length){full_window_string.erase(full_window_string.begin());}


		//std::cout << full_window_string << std::endl;

		// FIND MINIMIZER K-MER

		std::string forward_window_minimizer_kmer = find_minimizer_kmer(full_window_string, window_number, minimizer_len, window_size);
		std::string reverse_window_minimizer_kmer = find_minimizer_kmer(reverse_complement_seqstr(full_window_string), window_number, minimizer_len, window_size);

		std::string canonical_window_minimizer_kmer = forward_window_minimizer_kmer;
		for (int ii = 0; ii < window_number*minimizer_len; ii+=1)
		{
			if (forward_window_minimizer_kmer.at(ii) == reverse_window_minimizer_kmer.at(ii)){continue;}
			if (forward_window_minimizer_kmer.at(ii) < reverse_window_minimizer_kmer.at(ii)){
				canonical_window_minimizer_kmer = forward_window_minimizer_kmer;
				break;
			}
			if (forward_window_minimizer_kmer.at(ii) > reverse_window_minimizer_kmer.at(ii)){
				canonical_window_minimizer_kmer = reverse_window_minimizer_kmer;
				break;
			}
		}

		if(previous_window_minimizer_kmer == canonical_window_minimizer_kmer){continue;}

		previous_window_minimizer_kmer = canonical_window_minimizer_kmer;

		if(kmer_map.count(canonical_window_minimizer_kmer) == 0){kmer_map[canonical_window_minimizer_kmer] = 0;}

		kmer_map[canonical_window_minimizer_kmer] += 1;


		windows_checked+=1;
		if (windows_checked % 10000 == 0){std::cout << "Windows checked: " << windows_checked << std::endl;}
		
	}

	std::string write_file_path = output_file_path + ".txt";
	std::ofstream write_file(write_file_path);

	int found_kmer_count = 0;
	for(auto & pair : kmer_map)
	{
		found_kmer_count += 1;
		std::cout << pair.first << " " << pair.second << "\n";
		//std::cout << "k-mer: " << pair.first << " count: " << pair.second << "\n";
		write_file << pair.first << " " << pair.second << "\n";



	}
	std::cout << "k-mers found: " << found_kmer_count << std::endl;


	//std::cout << "Program run finished successfully. Thank you for using LoMeX." << std::endl;
	return 0;

	*/


}











