#include <iostream>
#include <fstream>

#include "spoa/spoa.hpp"

// To compile: g++ example2.cpp -std=c++11 -Iinclude/ -Lbuild/lib/ -lspoa -o example2 -g

// Read sequences from a fasta file and construct POA for them
// Arguments: fasta file
int main(int argc, char** argv) {

  std::vector<std::string> sequences;
  std::ifstream infile(argv[1]);
  std::string line;
  while(std::getline(infile, line)) {
    if (line[0] != '>')
      sequences.emplace_back(line);
  }
  
  auto alignment_engine = spoa::AlignmentEngine::Create(
    spoa::AlignmentType::kNW, 3, -5, -3);  // linear gaps
  //spoa::AlignmentType::kNW, 4, -2, -12,-6);  // affine gaps

  spoa::Graph graph{};

  for (const auto& it : sequences) {
    auto alignment = alignment_engine->Align(it, graph);
    graph.AddAlignment(alignment, it);
  }

  // Generate a single consensus sequence
  //auto consensus = graph.GenerateConsensus();

  //std::cerr << ">Consensus LN:i:" << consensus.size() << std::endl
  //          << consensus << std::endl;

  // Generate a set of consensus sequences
  std::vector<std::string> cons;
  int *bundles = new int[sequences.size()];
  graph.GenerateConsensi(&cons, bundles, 0.8);

  // Add the consensus sequences to POA
  for(auto c = cons.begin(); c != cons.end(); ++c) {
    std::cerr << "Consensus: " << (*c) << std::endl;
    auto alignment = alignment_engine->Align((*c), graph);
    graph.AddAlignment(alignment, (*c));
  }

  // Generate a multile alignment
  auto msa = graph.GenerateMultipleSequenceAlignment();

  int j = 0;
  for (const auto& it : msa) {
    if (j < sequences.size())
      std::cerr << bundles[j++] << "  " << it << std::endl;
    else
      std::cerr << (j++-sequences.size()+1) << "* " << it << std::endl;
      
  }

  return 0;
}
