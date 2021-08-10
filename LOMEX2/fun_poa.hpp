#ifndef FUN_POA_H
#define FUN_POA_H

#include <string>
#include <vector>


extern "C" 
{
 	#include "lpo.h"
	#include "msa_format.h"
	#include "align_score.h"
}



class LOMEXPOAGlue
{
private:
	ResidueScoreMatrix_T score_matrix;
	int score_matrix_int;
	std::vector<std::string> consensus_sequences;

public:

	void initialize_me(std::string score_matrix_path);

	std::vector<std::string>& get_consensus_sequences(){return consensus_sequences;}

	void poa_consensus_sequences_unrefined(std::vector<std::string>& input_sequences, float bundling_threshold, float support_threshold, int min_bundle_size_threshold, int min_char_support_threshold);

	int poa_consensus_sequences(Sequence_T **mainseqs, std::vector<std::string>& input_sequences, float bundling_threshold, float support_threshold, int min_bundle_size_threshold, int min_char_support_threshold);

};


void do_poa_stuff();

//int construct_seqs(Sequence_T **mainseqs, std::vector<std::string>& input_sequences, int do_switch_case);

std::tuple<int,uint64_t,uint64_t, uint64_t, uint64_t, uint64_t, uint64_t, uint64_t> construct_consensus_seqs(std::vector<std::string>& occurrence_sequences, std::vector<std::string>& consensus_sequences, ResidueScoreMatrix_T *score_matrix, int do_switch_case,
	float bundling_threshold, float support_threshold, int min_bundle_size_threshold, int min_char_support_threshold);

//std::vector<std::string> poa_consensus_sequences(std::vector<std::string> input_sequences, std::string score_matrix_path, float bundling_threshold);

#endif
