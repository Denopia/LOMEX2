
#include <iostream>
#include <tuple>
#include "kmer_index.hpp"
#include "fun_kmers.hpp"





void KMerIndex::initialize_variables(std::string input_kmers_path, int input_kmer_len)
{
	kmer_len = input_kmer_len;
	kmers_n = 0;
	occurrences_n = 0;
	kmers_path = input_kmers_path;
}


void KMerIndex::initialize_index(int minimizer_set_occ_threshold)
{
	std::ifstream kmer_file(kmers_path);
	std::string one_line;
	uint64_t kmerint;
	
	std::string minimizer_set;
	int minimizer_set_occ;

	int too_smalls = 0;
	int big_enoughs = 0;

	std::cout << "Doing minimizer set reading\n";

	while (std::getline(kmer_file, one_line))
	{
		if ((too_smalls+big_enoughs)%10000 == 0){std::cout << "k-mers read: " << too_smalls+big_enoughs <<  "\n";}
		
		tie(minimizer_set, minimizer_set_occ) = interpret_minimizer_set_line(one_line);

		if (minimizer_set_occ < minimizer_set_occ_threshold){
			too_smalls += 1;
			continue;
		} else {
			big_enoughs += 1;
			kmerint = map_str2int_small(minimizer_set);	
			relevant_kmers.insert(kmerint);
		}
	}
	std::cout << "Too smalls: " << too_smalls << "\n";
	std::cout << "Big enoughs: " << big_enoughs << "\n";
}

void KMerIndex::add_occurrence(uint64_t kmer, int read_id, int position, int length, int strand)
{
	if (relevant_kmers.count(kmer)>0)
	{
		kmers2occurrences[kmer].push_back(read_id);
		kmers2occurrences[kmer].push_back(position);
		kmers2occurrences[kmer].push_back(length);
		kmers2occurrences[kmer].push_back(strand);
	}
	else
	{
		std::cout << "## ERROR ## The k-mer you are trying to insert is not relevant\n";
	}
}

bool KMerIndex::is_this_relevant_kmer(uint64_t kmer)
{
	if (relevant_kmers.count(kmer)>0)
	{
		return true;
	}
	else
	{
		return false;
	}
}


uint64_t KMerIndex::canonical_kmer(uint64_t kmer)
{
	return kmer;

}

void KMerIndex::print_index()
{

}


// Basically a split function but probably a super slow one
std::tuple<std::string, int> interpret_minimizer_set_line(std::string one_line)
{
	std::string minimizer_set_string = "";
	std::string minimizer_set_occ_string = "";
	int part = 0;

	for (char c : one_line)
	{
		if (!std::isspace(c))
		{
			if (part == 0){minimizer_set_string = minimizer_set_string + c;}
			else if (part == 1){minimizer_set_occ_string = minimizer_set_occ_string + c;}
			else {break;}
		} else {
			part+=1;
		}
	}

	//std::cout << "Minimizer set: " << minimizer_set_string << "\n";
	//std::cout << "Occurrences: " << minimizer_set_occ_string << "\n";


	return std::make_tuple(minimizer_set_string, stoi(minimizer_set_occ_string));
}


