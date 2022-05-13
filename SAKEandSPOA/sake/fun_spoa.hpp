#ifndef FUN_SPOA_H
#define FUN_SPOA_H

#include <string>
#include <vector>


void construct_consensus_seqs(std::vector<std::string>& occurrence_sequences, std::vector<std::string>& consensus_sequences, 
								float bundling_threshold, float support_threshold, int min_bundle_size_threshold, int min_char_support_threshold);


#endif
