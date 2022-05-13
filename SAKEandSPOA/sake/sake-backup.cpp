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
#include "fun_spoa.hpp"
#include <omp.h>


namespace po = boost::program_options;
namespace bfs = boost::filesystem;
//using namespace std;
//using namespace boost::filesystem;


int main(int argc, char *argv[])
{

	/*
		Initialize variables for user given arguments.
	*/
	std::string reads_path, kmers_path, output_path, score_matrix_path;
	int window_number, window_size, window_offset, minimizer_len;
	float placeholder_float;

	/*
		Parse arguments given by the user.
	*/
	po::options_description desc("LOMEXII options");
	desc.add_options()
		("help,h", "Give help")
		("kmers,k", po::value<std::string>(& kmers_path)->default_value("na0"), "Path to the k-mers file")
		("reads,r", po::value<std::string>(& reads_path)->default_value("na1"), "Path to the reads file")
		("scores,m", po::value<std::string>(& score_matrix_path)->default_value("na3"), "Path to the score matrix file")
		("output,o", po::value<std::string>(& output_path)->default_value("na2"), "Path to the output file")
		("windows,n", po::value<int>(& window_number)->default_value(1), "Number of windows")
		("offset,f", po::value<int>(& window_offset)->default_value(1), "Window offset")
		("window-size,s", po::value<int>(& window_size)->default_value(1), "Size of window")
		("minimizer-len,l", po::value<int>(& minimizer_len)->default_value(1), "Minimizer length")
		("pl-float,q", po::value<float>(& placeholder_float)->default_value(1.1), "Placeholder float");
		
	po::variables_map vm;
	po::store(po::parse_command_line(argc, argv, desc), vm);
	po::notify(vm);


	// SOME HARD CODED PARAMETERS THAT NEED TO BE MADED MODIFIABLE WHEN RUNNING THE PROGRAM
	float bundling_threshold = 0.8;
	float support_threshold = 0.4;
	int minimizer_set_occ_threshold = 3;
	int min_bundle_size_threshold = 3;
	int min_char_support_threshold = 3;

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
	
	int minimizer_kmer_len = window_number * minimizer_len;

	KMerIndex kmer_index;
	kmer_index.initialize_variables(kmers_path, minimizer_kmer_len);
	kmer_index.initialize_index(minimizer_set_occ_threshold);


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
	
	std::vector<uint64_t> current_read_kmers_as_ints; //all minimizer length k-mers as integers, forward
	std::vector<uint64_t> current_read_kmers_as_ints_r; //all minimizer length k-mers as integers, reversed

	std::vector<uint64_t> current_read_minimizers; // all minimizers
	std::vector<uint64_t> current_read_minimizers_r;
	
	std::vector<uint64_t> current_read_minimizer_kmers; // all minimizers
	std::vector<uint64_t> current_read_minimizer_kmers_r;

	//std::vector<uint64_t> current_read_kmer_hashes;
	//std::vector<uint64_t> current_read_kmer_hashes_r;


	std::vector<int> current_read_minimizer_kmer_strands; 
		
	//std::vector<int> current_read_minimizer_position_global;
	//std::vector<int> current_read_minimizer_position_r_global;

	std::vector<int> current_read_minimizer_position_in_window;
	std::vector<int> current_read_minimizer_position_in_window_r;

	std::vector<std::vector<int>> current_read_minimizer_kmer_position_in_window_coordinates;
	//std::vector<std::vector<int>> current_read_minimizer_kmer_position_global_coordinates;
	

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
	std::cout << "Start reading reads\n";

	while(reads_file_reader.get_reads_left())
	{
		one_read = "";//-------------------------------------------------------------------------------------------------------------------------------------
		// === Read the next read ===
		reads_file_reader.roll_to_next_read();
		// === Clear variables related to the previous read ===
		current_read_kmers_as_ints.clear();
		current_read_kmers_as_ints_r.clear();

		current_read_minimizers.clear();
		current_read_minimizer_kmers.clear();
		current_read_minimizer_kmer_strands.clear();

		current_read_minimizers_r.clear();
		current_read_minimizer_kmers_r.clear();

		//current_read_kmer_hashes.clear();
		//current_read_kmer_hashes_r.clear();

		//current_read_minimizer_position_global.clear();
		//current_read_minimizer_position_r_global.clear();

		current_read_minimizer_position_in_window.clear();
		current_read_minimizer_position_in_window_r.clear();

		current_read_minimizer_kmer_position_in_window_coordinates.clear();
		//current_read_minimizer_kmer_position_global_coordinates.clear();

		
		current_kmer_as_int = 0;
		characters_read_from_read = 0;
		current_kmer_as_string = "";

		if (reads_file_reader.get_current_read_number() % 100 == 0){
			std::cout << "This is read number " << reads_file_reader.get_current_read_number() << "/" << reads_file_reader.get_total_reads() << " with length " << reads_file_reader.get_read_length() << "\n";
			//break;
		}
		//if (reads_file_reader.get_current_read_number() < 24000){continue;}

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
				current_read_kmers_as_ints_r.push_back(reverse_complement_seqint_small(current_kmer_as_int, minimizer_len));
				//std::cout << "k-mer: " + current_kmer_as_string << " is " << current_kmer_as_int << "\n";
			}
		} // All k-mers in the current read have been found and stored as ints

		all_reads.push_back(one_read);

		/*
		uint64_t forward_hash;
		uint64_t reverse_hash;
		for (int i = 0; i < current_read_kmers_as_ints.size(); i++)
		{
			forward_hash = current_read_kmers_as_ints[i];
			reverse_hash = reverse_complement_seqint_small(current_read_kmers_as_ints[i], minimizer_len);
			current_read_kmer_hashes.push_back(forward_hash);
			current_read_kmer_hashes_r.push_back(reverse_hash);
		}

		std::cout << "REVERSE AS INTS LEN = " << current_read_kmers_as_ints_r.size() <<"\n";
		std::cout << "REVERSE HASHES  LEN = " << current_read_kmer_hashes_r.size() <<"\n";
		*/

		

		//std::cout << "k-mers found as ints successfully\n";
		//std::cout << "There are " << current_read_kmers_as_ints.size() << " k-mers in the read\n\n";

		// === After all k-mers are in the vector it is time to find all the minimizers ===
		// ================================================================================
		// ================================================================================

		// This is something that can be made faster if we do not check every window separately

		// ================================================================================
		// ================================================================================

		// === Set the minimizer in this variable ===
		uint64_t i_am_the_minimizer;
		int minimizer_position;																												// WHAT IS MINIMIZER POSITIONS SUPPOSED TO BE?
		// === Now go through all windows to find the minimizers ===
		for (int i = 0; i+window_size-1 < current_read_kmers_as_ints.size(); i++)
		{
			// === Set the minimizer to the maximum value so it will be updated (or at least tied) ===
			//i_am_the_minimizer = std::numeric_limits<uint64_t>::max();
			i_am_the_minimizer = current_read_kmers_as_ints[i];
			// === Set the minimizer position at the start of the window (it will be updated) ===
			minimizer_position = i;
			// === Now check every k-mer in the window ===
			int position_in_window = 1;
			int minimizer_position_in_window = 0;
			for (int j = i+1; j < i+window_size; j++)
			{
				// === If a new minimizer is found, update variables === 																									HERE WE CALL THE MINIMIZER HASH FUNCTION
				if (minimizer_hash(current_read_kmers_as_ints[j]) <= minimizer_hash(i_am_the_minimizer))
				{
					i_am_the_minimizer = current_read_kmers_as_ints[j];
					minimizer_position = j;
					minimizer_position_in_window = position_in_window;
				}
				position_in_window += 1;
			}
			// === After all k-mers are checked, push the found minimizer in the minimizer vector
			current_read_minimizers.push_back(i_am_the_minimizer);
			current_read_minimizer_position_in_window.push_back(minimizer_position_in_window);
			//current_read_minimizer_position_global.push_back(minimizer_position);
		}

		// DO THE SAME FOR REVERSE STRAND
		uint64_t i_am_the_minimizer_r;
		//int minimizer_position_r;
		// === Now go through all windows to find the minimizers ===
		for (int i = 0; i+window_size-1 < current_read_kmers_as_ints_r.size(); i++)
		{
			// === Set the minimizer to the maximum value so it will be updated (or at least tied) ===
			//i_am_the_minimizer_r = std::numeric_limits<uint64_t>::max();
			i_am_the_minimizer_r = current_read_kmers_as_ints_r[i];
			// === Set the minimizer position at the start of the window (it will be updated) ===
			//minimizer_position_r = i;
			// === Now check every k-mer in the window ===
			int position_in_window_r = 1;
			int minimizer_position_in_window_r = 0;
			for (int j = i+1; j < i+window_size; j++)
			{
				// === If a new minimizer is found, update variables ===
				if (minimizer_hash(current_read_kmers_as_ints_r[j]) <= minimizer_hash(i_am_the_minimizer_r))
				{
					i_am_the_minimizer_r = current_read_kmers_as_ints_r[j];
					//minimizer_position_r = j;
					minimizer_position_in_window_r = position_in_window_r;
				}
				position_in_window_r += 1;
			}
			// === After all k-mers are checked, push the found minimizer in the minimizer vector
			current_read_minimizers_r.push_back(i_am_the_minimizer_r);
			current_read_minimizer_position_in_window_r.push_back(minimizer_position_in_window_r);
			//current_read_minimizer_position_r_global.push_back(minimizer_position_r);
		}


		// === Now that all minimizers are found, we can find all minimizer-k-mers ===



		//std::cout << "Minimizers found successfully\n";
		//std::cout << "There are " << current_read_minimizers.size() << " minimizers in the vector\n\n";

		

		// === First clear the minimizer-k-mer vector of previous read ===
		
		//current_read_minimizer_kmers.clear();


		// =========================================================  IT IS LIKELY THERE IS A BUG HERE IN CREATING THE REVERSE K-MER!!!!!!!!!!!!!!!!!!!!!!!!

		// Put minimizer-k-mer here
		uint64_t forward_minimizer_kmer;
		uint64_t reverse_minimizer_kmer;
		// This is the distance between two minimizer windows
		//int minimizer_window_distance = minimizer_len + window_size - 1;

		std::vector<int> minimizer_kmer_coordinates_in_window;
		std::vector<int> minimizer_kmer_coordinates_in_window_r;

		//std::vector<int> minimizer_kmer_coordinates_global;
		//std::vector<int> minimizer_kmer_coordinates_r_global;


		std::vector<int> current_insert_minimizer_kmer_coordinates_in_window;
		//std::vector<int> current_insert_minimizer_kmer_coordinates_global;

		for (int i = 0; i + (window_number-1)*window_offset < current_read_minimizers.size(); i++)
		{

			//std::cout << "Minimizer at position: " << i << "\n";
			forward_minimizer_kmer = 0;
			reverse_minimizer_kmer = 0;

			minimizer_kmer_coordinates_in_window.clear();
			minimizer_kmer_coordinates_in_window_r.clear();

			//minimizer_kmer_coordinates_global.clear();
			//minimizer_kmer_coordinates_r_global.clear();
			

			for (int j = 0; j < window_number; j++)
			{
				forward_minimizer_kmer = forward_minimizer_kmer << 2*minimizer_len;
				forward_minimizer_kmer = forward_minimizer_kmer | current_read_minimizers[i+j*window_offset];
			
				
				reverse_minimizer_kmer = reverse_minimizer_kmer << 2*minimizer_len;
				reverse_minimizer_kmer = reverse_minimizer_kmer | current_read_minimizers_r[i+(window_number-1-j)*window_offset]; // NOT SURE IF RIGHT????????

				
				minimizer_kmer_coordinates_in_window.push_back(current_read_minimizer_position_in_window[i+j*window_offset]); // DIS DONE MODIFIED
				//minimizer_kmer_coordinates_global.push_back(current_read_minimizer_position_global[i+j*window_offset]); // DIS DONE MODIFIED

				//minimizer_kmer_coordinates_in_window_r.push_back(current_read_minimizer_position_in_window_r[i+(window_number-1-j)*window_offset]); // AND DIS
				minimizer_kmer_coordinates_in_window_r.push_back(current_read_minimizer_position_in_window_r[i+j*window_offset]);
				//minimizer_kmer_coordinates_r_global.push_back(current_read_minimizer_position_r_global[i+(window_number-1-j)*window_offset]); // AND DIS
			}

			if (forward_minimizer_kmer <= reverse_minimizer_kmer){
				current_read_minimizer_kmers.push_back(forward_minimizer_kmer);
				current_read_minimizer_kmer_strands.push_back(1);
				current_read_minimizer_kmer_position_in_window_coordinates.push_back(minimizer_kmer_coordinates_in_window);
				//current_read_minimizer_kmer_position_global_coordinates.push_back(minimizer_kmer_coordinates_global);

				//std::cout << "FORWARD\n";
			} else {
				current_read_minimizer_kmers.push_back(reverse_minimizer_kmer);
				current_read_minimizer_kmer_strands.push_back(-1);
				current_read_minimizer_kmer_position_in_window_coordinates.push_back(minimizer_kmer_coordinates_in_window_r);
				//current_read_minimizer_kmer_position_global_coordinates.push_back(minimizer_kmer_coordinates_r_global);
				//std::cout << "==============================================BACK\n";
			}
		}

		//std::cout << "minimizer-k-mers found successfully\n";
		//std::cout << "There are " << current_read_minimizer_kmers.size() << " minimizer-k-mers in the vector\n\n";


		// === After this we can finally put the minimizer-k-mer occurrences into the minimizer-k-mer index ===
		uint64_t previous_insert_minimizer_kmer = std::numeric_limits<uint64_t>::max();
		uint64_t current_insert_minimizer_kmer;
		int current_insert_minimizer_kmer_strand;
		current_insert_minimizer_kmer_coordinates_in_window.clear();
		//current_insert_minimizer_kmer_coordinates_global.clear();
		int iloop = 0;
		int jloop = 0;
		int same_minimizer_kmer_streak;
		//int minimizer_kmer_full_window_size = window_number*(minimizer_len + window_size - 1);
		//int full_window_length = (window_number - 1) * window_offset + (window_size + minimizer_len - 1);

		while(iloop < current_read_minimizer_kmers.size())
		{
			current_insert_minimizer_kmer = current_read_minimizer_kmers.at(iloop);
			current_insert_minimizer_kmer_strand = current_read_minimizer_kmer_strands.at(iloop);
			current_insert_minimizer_kmer_coordinates_in_window = current_read_minimizer_kmer_position_in_window_coordinates.at(iloop);
			//current_insert_minimizer_kmer_coordinates_global = current_read_minimizer_kmer_position_global_coordinates.at(iloop);

			same_minimizer_kmer_streak = 0;
			jloop = iloop+1;
			while (jloop < current_read_minimizer_kmers.size())
			{
				if (current_insert_minimizer_kmer != current_read_minimizer_kmers.at(jloop)){break;}
				if (current_insert_minimizer_kmer_strand != current_read_minimizer_kmer_strands.at(jloop)){break;}                                           // IS THIS NECESSARY?
				//std::cout << "Extended minimizer-k-mer\n";
				jloop+=1;
				same_minimizer_kmer_streak+=1;
			}

			if (kmer_index.is_this_relevant_kmer(current_insert_minimizer_kmer))
			{
				//std::cout << "Inserted minimizer-k-mer\n";
				//kmer_index.add_occurrence(current_insert_minimizer_kmer, reads_file_reader.get_current_read_number(), 
				//							iloop, minimizer_kmer_full_window_size + same_minimizer_kmer_streak, 
				//							current_insert_minimizer_kmer_strand, current_insert_minimizer_kmer_coordinates_in_window);


				// Start position is the current position plus the window position of the first minimizer
				int input_start_position_for_poa = iloop + current_insert_minimizer_kmer_coordinates_in_window[0];
				// Global position of the start of the last minimizer
				int position_of_last_minimizer = iloop + (window_number-1)*window_offset + current_insert_minimizer_kmer_coordinates_in_window[window_number-1];
				// Global end position of the last minimizer
				int end_of_last_minimizer = position_of_last_minimizer + minimizer_len - 1;
				// Length of the area covered by this minimizer-k-mer (so that tails are first and last minimizer)
				int input_length_for_spoa = end_of_last_minimizer - input_start_position_for_poa + 1;

				kmer_index.add_occurrence(current_insert_minimizer_kmer, reads_file_reader.get_current_read_number(), 
											input_start_position_for_poa, input_length_for_spoa, current_insert_minimizer_kmer_strand);

				//kmer_index.add_occurrence(current_insert_minimizer_kmer, reads_file_reader.get_current_read_number(), 
				//							iloop, full_window_length + same_minimizer_kmer_streak, 
				//							current_insert_minimizer_kmer_strand);
			} 
			//else {
			//	std::cout << "IRRELEVANT MINIMIZER-K-MER ERROR (BFCOUNTER AND LOMEX DISAGREE)\n";
			//}
			iloop = jloop;
		}

		//std::cout << "\n";
	} // Every read has been read after this


	return 0;


	std::ofstream output_file;
	output_file.open(output_path);
	

	int number_of_kmers = kmer_index.get_kmers_n();

	//for (auto minikmer : kmer_index.get_kmer_set())
	// num_threads(8)
	omp_set_num_threads(1);
	#pragma omp parallel for
	for (int ompi = 0; ompi < number_of_kmers; ompi++)
	{
		if (ompi == 1){continue;}
		uint64_t minikmer = kmer_index.get_kmer_vector()[ompi];
		std::string minikmer_as_string = map_int2str_small(minikmer, window_number*minimizer_len); // NEW 
		int minimizer_kmer_id = ompi;

		//std::cout << "Now starting to work with minimizer-k-mer ID (or rather its number):" << minimizer_kmer_id << " as string: " << minikmer_as_string << "\n";
		
		std::vector<std::string> occurrence_strings;
		std::vector<std::string> consensus_sequences;

		//std::cout << "=========================================================\n";
		//std::cout << "K-MER " << minimizer_kmer_id << " " << map_int2str_small(minikmer, minimizer_kmer_len) << " OCCURRENCES\n";
		int idk = 0;
		
		// Create sequence anchor map
		int nof_seqs = kmer_index.get_nof_occurrences(minikmer);

		if (nof_seqs == 0){continue;}

		int instances = 0;
		
		while (idk < kmer_index.get_kmer_info(minikmer).size())
		{
			//std::cout << "IDK START = " << idk << "\n";
			int krid = kmer_index.get_kmer_info(minikmer)[idk];
			idk+=1;
			int krpos = kmer_index.get_kmer_info(minikmer)[idk];
			idk+=1;
			int krlen = kmer_index.get_kmer_info(minikmer)[idk];
			idk+=1;
			int krstrand = kmer_index.get_kmer_info(minikmer)[idk];
			idk+=1;

			//std::cout << "IDK END = " << idk << "\n";
			//std::cout << "READ ID: " << krid << "\n";
			//std::cout << "POSITION IN READ: " << krpos << "\n";
			//std::cout << "LENGTH: " << krlen << "\n";
			//std::cout << "STRAND: " << krstrand << "\n";
			//std::cout << "AS STRING: " << all_reads[krid-1].substr(krpos, krlen) << "\n";
			//std::cout << "ANCHOR POSITIONS IN WINDOWS: ";	
			
			// PUT THE CORRECT ORIENTATION OF THE OCCURRENCE TO THE ARRAY
			if (krstrand == 1){
				occurrence_strings.push_back(all_reads[krid-1].substr(krpos, krlen));
			} else {
				occurrence_strings.push_back(reverse_complement_seqstr(all_reads[krid-1].substr(krpos, krlen)));
			}
			instances = instances + 1;
		}

		//std::cout << "This minimizer-k-mer has " << instances << " instances in the reads which are: \n";
		//for (int pp = 0; pp < occurrence_strings.size(); pp++){
		//	std::cout << occurrence_strings[pp] << "\n";
		//}
		//std::cout << "those were the strings\n\n";

		if (occurrence_strings.size() == 0){
			occurrence_strings.clear();
			consensus_sequences.clear();
			continue;
		}

		construct_consensus_seqs(occurrence_strings, consensus_sequences, bundling_threshold, support_threshold, min_bundle_size_threshold, min_char_support_threshold);

		//std::cout << "Construct done\n";
		//std::cout << "CONSENSUS SEQUENCES FOR " <<  map_int2str_small(minikmer, minimizer_kmer_len) << " :\n";
		
		#pragma omp critical(print_output)
		{
			std::cout << "This was minimizer-k-mer ID:" << minimizer_kmer_id << " as string: " << minikmer_as_string << "\n";
			for (auto css : consensus_sequences)
			{
				//std::cout << css << " " << minimizer_kmer_id << "\n";
				output_file << css << " " << minimizer_kmer_id << " " << minikmer_as_string << "\n";
			}
			std::cout << "----------------------------------------------------------------------------------------------\n";
		}

		occurrence_strings.clear();
		consensus_sequences.clear();


	

	}
	
	output_file.close();

	return 0;

}











