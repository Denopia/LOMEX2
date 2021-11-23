#include <iostream>

#include "spoa/spoa.hpp"

// To compile: g++ example.cpp -std=c++11 -Iinclude/ -Lbuild/lib/ -lspoa -o example -g

// Measure the runtime of generating a POA of a set of sequences
int main(int argc, char** argv) {

  // The sequences
  std::vector<std::string> sequences = {
      "CATAAAAGAACGTAGGTCGCCCGTCCGTAACCTGTCGGATCACCGGAAAGGACCCGTAAAGTGATAATGAT",
      "ATAAAGGCAGTCGCTCTGTAAGCTGTCGATTCACCGGAAAGATGGCGTTACCACGTAAAGTGATAATGATTAT",
      "ATCAAAGAACGTGTAGCCTGTCCGTAATCTAGCGCATTTCACACGAGACCCGCGTAATGGG",
      "CGTAAATAGGTAATGATTATCATTACATATCACAACTAGGGCCGTATTAATCATGATATCATCA",
      "GTCGCTAGAGGCATCGTGAGTCGCTTCCGTACCGCAAGGATGACGAGTCACTTAAAGTGATAAT",
      "CCGTAACCTTCATCGGATCACCGGAAAGGACCCGTAAATAGACCTGATTATCATCTACAT"
  };

  // Generate the same POA 100000 times
  for(int i = 0; i < 100000; i++) {
    
    auto alignment_engine = spoa::AlignmentEngine::Create(
      spoa::AlignmentType::kNW, 3, -5, -3);  // linear gaps

    spoa::Graph graph{};

    // Add the sequences to POA
    for (const auto& it : sequences) {
      auto alignment = alignment_engine->Align(it, graph);
      graph.AddAlignment(alignment, it);
    }

    // Generate consensus sequences
    std::vector<std::string> cons;
    int bundles[6];
    graph.GenerateConsensi(&cons, bundles, 0.75);

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

    //int j = 0;
    //for (const auto& it : msa) {
    // std::cerr << bundles[j++] << " " << it << std::endl;
    //}
  }

  return 0;
}
