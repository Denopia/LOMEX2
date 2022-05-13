#include <iostream>
#include <vector>
#include <stdio.h>
#include <string>
#include <typeinfo>
#include <unordered_map>
#include "fun_spoa.hpp"
#include "spoa.hpp"
#include "fun_kmers.hpp"





std::string randomize_ns(std::string s){
  std::string s2 = "";
  for (char c : s){
    if (c == 'N'){
      s2 = s2 + random_nuc();
    } else {
      s2 = s2 + c;
    }
  }
  return s2;
}


// THE MAIN FUNCTION THAT IS USED TO CONSTRUCT CONSENSUS SEQUENCES
void construct_consensus_seqs(std::vector<std::string>& occurrence_sequences, std::vector<std::string>& consensus_sequences, 
															float bundling_threshold, float support_threshold, int min_bundle_size_threshold, int min_char_support_threshold)
{

  //std::cout << "STARTING TO BUILD ALIGNMENT ENGINE\n";

  auto alignment_engine = spoa::AlignmentEngine::Create(spoa::AlignmentType::kNW, 3, -5, -3);  // linear gaps
  // match, mismatch, gap (linear, open and extend sama), fourth parameter affine= gap extend



  spoa::Graph graph{};

  //std::cout << "PUTTING STRINGS IN ALIGNMENT\n";

  // Add the sequences to POA
  for (const auto& it : occurrence_sequences) {
    std::string itr = randomize_ns(it);
    auto alignment = alignment_engine->Align(it, graph);
    graph.AddAlignment(alignment, it);
  }

  // Generate consensus sequences
  std::vector<std::string> cons;
  //int bs = 6;
  int bundles[occurrence_sequences.size()]; // IS THIS CORRECT?
  //int bundling_threshold = 0.75

  //std::cout << "START BUILDING CONSENSUS\n";

  graph.GenerateConsensi(&cons, bundles, bundling_threshold);

  //std::cout << "CONSENSUS CREATED\n";

  //for(auto c = cons.begin(); c != cons.end(); ++c) {
  //  std::cerr << "Consensus: " << (*c) << std::endl;
  //}

  //auto consensus = graph.GenerateConsensus();
  
  //std::cerr << ">Consensus LN:i:" << consensus.size() << std::endl
  //          << consensus << std::endl;

  // Add the consensus sequences to POA
  for(auto c = cons.begin(); c != cons.end(); ++c) {
    auto alignment = alignment_engine->Align((*c), graph);
    graph.AddAlignment(alignment, (*c));
  }

  // Generate MSA
  auto msa = graph.GenerateMultipleSequenceAlignment();

  //std::cout << "MSA CREATED\n";

  std::vector<std::string> all_msa_seqs;

  int j = 0;
  for (const auto& it : msa) {
    all_msa_seqs.push_back(it);
    /*
    if (j < occurrence_sequences.size()){
      std::cout << bundles[j++] << " " << it << std::endl;
    } else {
      std::cout << "X " << it << std::endl;
      consensus_sequences.push_back(it);
    }
    */
  }




  std::unordered_map<int,int> bundle_counts;
  std::vector<std::string> input_seqs_aligned;
  std::unordered_map<int, std::vector<std::string>> input_seq_bundles;
  std::vector<std::string> consensus_seqs_aligned;
  std::string conseq;

  int concount = 1;
  int seq_bundle;
  
  //std::cout << "FILTERING STARTS\n";

  for (int i = 0; i < all_msa_seqs.size(); i += 1)
  {
    conseq = all_msa_seqs[i];
    //std::cout << "*****\n";
    //std::cout << "SEQ NAME: " << lpo_out->source_seq[i].name << "\n";
    //std::cout << "SEQ TITLE: " << lpo_out->source_seq[i].title << "\n";
    //std::cout << "SEQ BUNDLE: " << lpo_out->source_seq[i].bundle_id << "\n";
    
    // Get sequence's bundle
    if (i < occurrence_sequences.size()){
      seq_bundle = bundles[i];
    } else {
      seq_bundle = concount;
      concount += 1;
    }
    
    // Input sequence bundles
    if (i < occurrence_sequences.size())
    {
      if (seq_bundle > 0)
      {
        if (bundle_counts.count(seq_bundle) == 0){
          bundle_counts[seq_bundle] = 0;
        }
        bundle_counts[seq_bundle] += 1;
        input_seq_bundles[seq_bundle].push_back(conseq);
      }   
      //continue;
    }

    // Print debug stuff
    /*
    if (i < occurrence_sequences.size()){
      std::cout << seq_bundle << "\t" << conseq << std::endl;
    } 
    */

    // Determine if consensus sequence is good enough to be put into the 
    // This must be a consensus sequence

    if (i >= occurrence_sequences.size())
    {
      bool is_real_bundle = bundle_counts.count(seq_bundle) > 0;
      bool is_supported_percent = bundle_counts[seq_bundle] > support_threshold * occurrence_sequences.size();
      bool is_supported_hard = bundle_counts[seq_bundle] > min_bundle_size_threshold;

      bool enough_support = true;
      std::vector<int> character_supports;

      for (int n = 0; n < all_msa_seqs[i].length(); n+=1)
      {
        int position_support = 0;

        // Deletions (dashses) do not need to be supported
        if (conseq.at(n) == '-')
        {
          character_supports.push_back(0);
          continue;
        }

        for (int m = 0; m < input_seq_bundles[seq_bundle].size(); m+=1)
        {
          if (input_seq_bundles[seq_bundle][m].at(n) == conseq.at(n)){position_support += 1;}
        }

        if (position_support < min_char_support_threshold)
        {
          character_supports.push_back(-1);
          //enough_support = false;
          //break;
        } else {
          character_supports.push_back(1);
        }
      }


      // Find supported start
      int supported_start = 0;
      for (int csi = 0; csi < conseq.length(); csi+=1)
      {
        if (character_supports.at(csi) == 1){supported_start = csi;break;}
      }

      // Find supported end
      int supported_end = 0;
      for (int csi = conseq.length()-1; csi >= 0; csi-=1)
      {
        if (character_supports.at(csi) == 1){supported_end = csi;break;}
      }
      
      // No random unsupported characters in the middle
      int from_support_to_no_support = 0;
      int from_no_support_to_support = 0;
      int illegal_dip_in_no_support = 0;
      int last_non_deletion = character_supports.at(0);
      for (int csi = 1; csi < conseq.length(); csi+=1)
      {
        if ((character_supports.at(csi) == -1) && (last_non_deletion == 1)){from_support_to_no_support+=1;}
        if ((character_supports.at(csi) == 1) && (last_non_deletion == -1)){
          from_no_support_to_support+=1;
          if (from_support_to_no_support > 0){illegal_dip_in_no_support += 1;}
        }
        if (character_supports.at(csi) == -1){last_non_deletion = -1;}
        if (character_supports.at(csi) == 1){last_non_deletion = 1;}
      }

      // Determine if enough support
      if (from_support_to_no_support > 1){enough_support = false;}
      if (from_no_support_to_support > 1){enough_support = false;}
      if (illegal_dip_in_no_support > 0){enough_support = false;}

      // Print accordingly
      if (enough_support && is_real_bundle && is_supported_percent && is_supported_hard){

        //consensus_sequences.push_back(conseq);
        consensus_sequences.push_back(conseq.substr(supported_start, supported_end-supported_start+1));
        //std::cout << seq_bundle << "\t" << conseq << " ** ("<< supported_start << "," << supported_end << ")\n";
        //std::cout << lpo_out->source_seq[i].bundle_id << "  " << conseq << " ** ("<< supported_start << "," << supported_end << ")\n";
      } else {
        ;
        //std::cout << seq_bundle << "\t" << conseq << " *\n" ;
        //std::cout << lpo_out->source_seq[i].bundle_id << "  " << conseq << " * \n";
      }
      character_supports.clear();
    }
  }


}
