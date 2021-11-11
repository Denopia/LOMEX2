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


std::tuple<int,uint64_t,uint64_t, uint64_t, uint64_t, uint64_t, uint64_t, uint64_t> construct_consensus_seqs(int split_align, int aggro, int prog, int glob, int nof_anchors, int **anchors, int minkmerid, std::vector<std::string>& occurrence_sequences, std::vector<std::string>& consensus_sequences, 
																												ResidueScoreMatrix_T *score_matrix, int do_switch_case, float bundling_threshold, float support_threshold, 
																												int min_bundle_size_threshold, int min_char_support_threshold);



#endif
