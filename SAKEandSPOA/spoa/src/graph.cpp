// Copyright (c) 2020 Robert Vaser

#include "spoa/graph.hpp"

#include <algorithm>
#include <cassert>
#include <fstream>
#include <stack>
#include <stdexcept>
#include <unordered_set>

namespace spoa {

Graph::Node::Node(std::uint32_t id, std::uint32_t code)
    : id(id),
      code(code),
      inedges(),
      outedges(),
      aligned_nodes() {
}

Graph::Node* Graph::Node::Successor(std::uint32_t label) const {
  for (const auto& it : outedges) {
    auto jt = std::find(it->labels.begin(), it->labels.end(), label);
    if (jt != it->labels.end()) {
      return it->head;
    }
  }
  return nullptr;
}

std::uint32_t Graph::Node::Coverage() const {
  std::unordered_set<std::uint32_t> labels;
  for (const auto& it : inedges) {
    std::copy(
        it->labels.begin(),
        it->labels.end(),
        std::inserter(labels, labels.end()));
  }
  for (const auto& it : outedges) {
    std::copy(
        it->labels.begin(),
        it->labels.end(),
        std::inserter(labels, labels.end()));
  }
  return labels.size();
}

Graph::Edge::Edge(
    Node* tail,
    Node* head,
    std::uint32_t label,
    std::uint32_t weight)
    : tail(tail),
      head(head),
      labels(1, label),
      weight(weight) {
}

void Graph::Edge::AddSequence(std::uint32_t label, std::uint32_t w) {
  labels.emplace_back(label);
  weight += w;
}

Graph::Graph()
    : num_codes_(0),
      coder_(256, -1),
      decoder_(256, -1),
      sequences_(),
      nodes_(),
      edges_(),
      rank_to_node_(),
      consensus_() {
}

Graph::Node* Graph::AddNode(std::uint32_t code) {
  nodes_.emplace_back(new Node(nodes_.size(), code));
  return nodes_.back().get();
}

void Graph::AddEdge(Node* tail, Node* head, std::uint32_t weight) {
  for (const auto& it : tail->outedges) {
    if (it->head == head) {
      it->AddSequence(sequences_.size(), weight);
      return;
    }
  }
  edges_.emplace_back(new Edge(tail, head, sequences_.size(), weight));
  tail->outedges.emplace_back(edges_.back().get());
  head->inedges.emplace_back(edges_.back().get());
}

Graph::Node* Graph::AddSequence(
    const char* sequence,
    const std::vector<std::uint32_t>& weights,
    std::uint32_t begin,
    std::uint32_t end) {
  if (begin == end) {
    return nullptr;
  }
  Node* prev = nullptr;
  for (std::uint32_t i = begin; i < end; ++i) {
    auto curr = AddNode(coder_[sequence[i]]);
    if (prev) {  // both nodes contribute to the weight
      AddEdge(prev, curr, weights[i - 1] + weights[i]);
    }
    prev = curr;
  }
  return nodes_[nodes_.size() - (end - begin)].get();
}

void Graph::AddAlignment(
    const Alignment& alignment,
    const std::string& sequence,
    std::uint32_t weight) {
  AddAlignment(alignment, sequence.c_str(), sequence.size(), weight);
}

void Graph::AddAlignment(
    const Alignment& alignment,
    const char* sequence, std::uint32_t sequence_len,
    std::uint32_t weight) {
  std::vector<std::uint32_t> weights(sequence_len, weight);
  AddAlignment(alignment, sequence, sequence_len, weights);
}

void Graph::AddAlignment(
    const Alignment& alignment,
    const std::string& sequence,
    const std::string& quality) {
  AddAlignment(
      alignment,
      sequence.c_str(), sequence.size(),
      quality.c_str(), quality.size());
}

void Graph::AddAlignment(
    const Alignment& alignment,
    const char* sequence, std::uint32_t sequence_len,
    const char* quality, std::uint32_t quality_len) {
  std::vector<std::uint32_t> weights;
  for (std::uint32_t i = 0; i < quality_len; ++i) {
    weights.emplace_back(quality[i] - 33);  // Phred quality
  }
  AddAlignment(alignment, sequence, sequence_len, weights);
}

void Graph::AddAlignment(
    const Alignment& alignment,
    const std::string& sequence,
    const std::vector<std::uint32_t>& weights) {
  AddAlignment(alignment, sequence.c_str(), sequence.size(), weights);
}

void Graph::AddAlignment(
    const Alignment& alignment,
    const char* sequence, std::uint32_t sequence_len,
    const std::vector<std::uint32_t>& weights) {
  if (sequence_len == 0) {
    return;
  }
  if (sequence_len != weights.size()) {
    throw std::invalid_argument(
        "[spoa::Graph::AddAlignment] error: "
        "sequence and weights are of unequal size!");
  }

  for (std::uint32_t i = 0; i < sequence_len; ++i) {
    if (coder_[sequence[i]] == -1) {
      coder_[sequence[i]] = num_codes_;
      decoder_[num_codes_++] = sequence[i];
    }
  }

  if (alignment.empty()) {
    sequences_.emplace_back(AddSequence(sequence, weights, 0, sequence_len));
    TopologicalSort();
    return;
  }

  std::vector<std::uint32_t> valid;
  for (const auto& it : alignment) {
    if (it.second != -1) {
      if (it.second < 0 || it.second >= static_cast<std::int32_t>(sequence_len)) {  // NOLINT
        throw std::invalid_argument(
            "[spoa::Graph::AddAlignment] error: invalid alignment");
      }
      valid.emplace_back(it.second);
    }
  }
  if (valid.empty()) {
    throw std::invalid_argument(
        "[spoa::Graph::AddAlignment] error: missing sequence in alignment");
  }

  // add unaligned bases
  Node* begin = AddSequence(sequence, weights, 0, valid.front());
  Node* prev = begin ? nodes_.back().get() : nullptr;
  Node* last = AddSequence(sequence, weights, valid.back() + 1, sequence_len);

  // add aligned bases
  for (const auto& it : alignment) {
    if (it.second == -1) {
      continue;
    }

    std::uint32_t code = coder_[sequence[it.second]];
    Node* curr = nullptr;
    if (it.first == -1) {
      curr = AddNode(code);
    } else {
      auto jt = nodes_[it.first].get();
      if (jt->code == code) {
        curr = jt;
      } else {
        for (const auto& kt : jt->aligned_nodes) {
          if (kt->code == code) {
            curr = kt;
            break;
          }
        }
        if (!curr) {
          curr = AddNode(code);
          for (const auto& kt : jt->aligned_nodes) {
            kt->aligned_nodes.emplace_back(curr);
            curr->aligned_nodes.emplace_back(kt);
          }
          jt->aligned_nodes.emplace_back(curr);
          curr->aligned_nodes.emplace_back(jt);
        }
      }
    }
    if (!begin) {
      begin = curr;
    }
    if (prev) {  // both nodes contribute to weight
      AddEdge(prev, curr, weights[it.second - 1] + weights[it.second]);
    }
    prev = curr;
  }
  if (last) {
    AddEdge(prev, last, weights[valid.back()] + weights[valid.back() + 1]);
  }
  sequences_.emplace_back(begin);

  TopologicalSort();
}

void Graph::TopologicalSort() {
  rank_to_node_.clear();

  std::vector<std::uint8_t> marks(nodes_.size(), 0);
  std::vector<bool> ignored(nodes_.size(), 0);

  std::stack<Node*> stack;
  for (const auto& it : nodes_) {
    if (marks[it->id] != 0) {
      continue;
    }
    stack.push(it.get());
    while (!stack.empty()) {
      auto curr = stack.top();
      bool is_valid = true;
      if (marks[curr->id] != 2) {
        for (const auto& jt : curr->inedges) {
          if (marks[jt->tail->id] != 2) {
            stack.push(jt->tail);
            is_valid = false;
          }
        }
        if (!ignored[curr->id]) {
          for (const auto& jt : curr->aligned_nodes) {
            if (marks[jt->id] != 2) {
              stack.push(jt);
              ignored[jt->id] = true;
              is_valid = false;
            }
          }
        }

        assert((is_valid || marks[curr->id] != 1) && "Graph is not a DAG");

        if (is_valid) {
          marks[curr->id] = 2;
          if (!ignored[curr->id]) {
            rank_to_node_.emplace_back(curr);
            for (const auto& jt : curr->aligned_nodes) {
              rank_to_node_.emplace_back(jt);
            }
          }
        } else {
          marks[curr->id] = 1;
        }
      }

      if (is_valid) {
        stack.pop();
      }
    }
  }

  assert(IsTopologicallySorted() && "Graph is not topologically sorted");
}

bool Graph::IsTopologicallySorted() const {
  assert(nodes_.size() == rank_to_node_.size() && "Topological sort not called ");  // NOLINT

  std::vector<bool> visited(nodes_.size(), 0);
  for (const auto& it : rank_to_node_) {
    for (const auto& jt : it->inedges) {
      if (!visited[jt->tail->id]) {
        return false;
      }
    }
    visited[it->id] = 1;
  }

  return true;
}

std::vector<std::uint32_t> Graph::InitializeMultipleSequenceAlignment(
    std::uint32_t* row_size) const {
  std::vector<std::uint32_t> dst(nodes_.size());
  std::uint32_t j = 0;
  for (std::uint32_t i = 0; i < rank_to_node_.size(); ++i, ++j) {
    auto it = rank_to_node_[i];
    dst[it->id] = j;
    for (const auto& jt : it->aligned_nodes) {
      dst[jt->id] = j;
      ++i;
    }
  }
  if (row_size) {
    *row_size = j;
  }
  return dst;
}

std::vector<std::string> Graph::GenerateMultipleSequenceAlignment(
    bool include_consensus) {
  std::uint32_t row_size = 0;
  auto node_id_to_column = InitializeMultipleSequenceAlignment(&row_size);

  std::vector<std::string> dst;
  for (std::uint32_t i = 0; i < sequences_.size(); ++i) {
    std::string row(row_size, '-');
    auto it = sequences_[i];
    while (true) {
      row[node_id_to_column[it->id]] = decoder_[it->code];
      if (!(it = it->Successor(i))) {
        break;
      }
    }
    dst.emplace_back(row);
  }
  if (include_consensus) {
    TraverseHeaviestBundle();
    std::string row(row_size, '-');
    for (const auto& it : consensus_) {
      row[node_id_to_column[it->id]] = decoder_[it->code];
    }
    dst.emplace_back(row);
  }

  return dst;
}

std::string Graph::GenerateConsensus() {
  TraverseHeaviestBundle();
  std::string dst{};
  for (const auto& it : consensus_) {
    dst += decoder_[it->code];
  }
  return dst;
}

std::string Graph::GenerateConsensus(
    std::vector<std::uint32_t>* summary,
    bool verbose) {
  if (!summary) {
    throw std::invalid_argument(
        "[spoa::Graph::GenerateConsensus] error: invalid ptr to summary");
  }

  auto dst = GenerateConsensus();

  summary->clear();
  if (!verbose) {
    for (const auto& it : consensus_) {
      summary->emplace_back(0);
      summary->back() += it->Coverage();
      for (const auto& jt : it->aligned_nodes) {
        summary->back() += jt->Coverage();
      }
    }
  } else {
    summary->resize((num_codes_ + 1) * consensus_.size(), 0);
    auto node_id_to_column = InitializeMultipleSequenceAlignment();

    for (std::uint32_t i = 0; i < sequences_.size(); ++i) {
      Node* it = sequences_[i];
      std::uint32_t c = 0, p, column = node_id_to_column[it->id];
      bool is_gap = false;
      while (true) {
        for (; c < consensus_.size(); ++c) {
          if (node_id_to_column[consensus_[c]->id] < column) {
            continue;
          } else {
            if (node_id_to_column[consensus_[c]->id] == column) {
              if (is_gap) {
                for (std::uint32_t j = p + 1; j < c; ++j) {
                  ++(*summary)[num_codes_ * consensus_.size() + j];
                }
              }
              is_gap = true;
              p = c;
              ++(*summary)[it->code * consensus_.size() + c];
            }
            break;
          }
        }
        if (c == consensus_.size() || !(it = it->Successor(i))) {
          break;
        }
        column = node_id_to_column[it->id];
      }
    }
  }

  return dst;
}

void Graph::GenerateConsensi(
    std::vector<std::string>* consensi,
    int* bundles, double bundling_threshold) {
  int absolute_minimum_good_edge_score = 3;
  
  int bundle_id = 1;
  int num_bundled = 0;
  for(int i = 0; i < sequences_.size(); i++) {
    bundles[i] = 0;
  }

  int *seq_len = new int[sequences_.size()];
  for(int i = 0; i < sequences_.size(); i++) {
    seq_len[i] = 0;
    for(auto n = sequences_[i]; n != nullptr; n = n->Successor(i)) {
      seq_len[i]++;
    }
  }
  
  while(num_bundled < sequences_.size()) {
    //std::cout << "Num bundled is: " << num_bundled << "\n";
    int avg_seq_len = 0;
    int num_seq = 0;
    for(int i = 0; i < sequences_.size(); i++) {
      if (!bundles[i]){
        num_seq++;
        avg_seq_len += seq_len[i];
      }
    }
    if (num_seq > 0) {
      avg_seq_len = avg_seq_len/num_seq;
    } else {
      avg_seq_len = 0;
    }

    if (avg_seq_len > edges_.size()) {
      avg_seq_len = edges_.size()-1;
    }

    //std::cout << "Avg: " << avg_seq_len << std::endl;
    
    // Classify edges into good and bad
    sort(edges_.begin(), edges_.end(), [](const std::unique_ptr<Edge> &lhs, const std::unique_ptr<Edge> &rhs)
    {
      return lhs->weight > rhs->weight;
    });


    if (avg_seq_len < 0) {avg_seq_len = 0;}
    if (avg_seq_len >= edges_.size()) {avg_seq_len = edges_.size() - 1;}

    //std::cout << "HERE\n";
    //std::cout << "EDGES SIZE " << edges_.size() << "\n";
    //std::cout << "AVG LEN " << avg_seq_len << "\n";
    int good_threshold = edges_[avg_seq_len]->weight;
    //std::cout << "AND THERE\n";
    if (good_threshold < absolute_minimum_good_edge_score){
      good_threshold = absolute_minimum_good_edge_score;
    }
    



    // Find the best sequence
    //std::cout << "In da loop\n";
    int num_good_edges = 0;
    for(auto it = edges_.begin(); it != edges_.end(); ++it) {
      //std::cout << "bloop\n";
      if ((*it)->weight >= good_threshold){
        num_good_edges++;
      }
    }
    int bestscore = 99999999;
    int best_seq = -1;
    //std::cerr << good_threshold << " " << num_good_edges << " " << edges_.size() << std::endl;
    //std::cerr << sequences_.size() << std::endl;
    for(int it = 0; it < sequences_.size(); it++) {
      if (bundles[it] != 0){
        continue;
      }
      int good_edges_in_seq = 0;
      Node *prev = sequences_[it];
      Node *next = nullptr;
      while(prev) {
        next = nullptr;
        for(const auto& it2 : prev->outedges) {
          auto jt = std::find(it2->labels.begin(), it2->labels.end(), it);
          if (jt != it2->labels.end()) {
            if (it2->weight >= good_threshold){
              good_edges_in_seq++;
            }
            next =it2->head;
            break;
          }
        }
        prev = next;
      }
      // Good edges not in this seq
      int current_score = num_good_edges - good_edges_in_seq;
      current_score = current_score + (seq_len[it]-1-good_edges_in_seq);
      //std::cerr << current_score << std::endl;
      if (current_score < bestscore){
        bestscore = current_score;
        best_seq = it;
      }
    }

    // If there is no best sequence, get out of here
    if(best_seq == -1){
      //std::cout << "No best seq\n";
      for(int ie = 0; ie < sequences_.size(); ie++) {
        if (bundles[ie] == 0){
          bundles[ie] = -1;
        }
      }
      delete [] seq_len;
      return;
    }
    //std::cerr << "Best seq: " << best_seq << std::endl;
  
    // Reweigh the best sequence
    int weight_adjust = 5;
    Node *prev = sequences_[best_seq];
    Node *next = nullptr;
    while(prev) {
      next = nullptr;
      for(const auto& it2 : prev->outedges) {
        auto jt = std::find(it2->labels.begin(), it2->labels.end(), best_seq);
        if (jt != it2->labels.end()) {
          it2->weight += weight_adjust;
          next =it2->head;
          break;
        }
      }
      prev = next;
    }

    // Find the bundle
    TraverseHeaviestBundle();
    std::vector<int> bundle;
    for(int it = 0; it < sequences_.size(); it++) {
      if (bundles[it] != 0){
        continue;
      }
      int current_score = 0;
      int good=0, all=0;
      for(auto n = sequences_[it]; n != nullptr; n = n->Successor(it)) {
        all++;
        auto ret = std::find(consensus_.begin(), consensus_.end(), n);
        if (ret != consensus_.end()) {
          good++;
        }
      }
      // std::cerr << good << " " << all << " " << (double)good/(double)all <<  std::endl;
      if ((double)good / (double)all >= bundling_threshold) {
        bundle.emplace_back(it);
      }
    }
  
    // Reweigh sequences in bundle
    for(auto it = bundle.begin(); it != bundle.end(); ++it) {
      if (*it != best_seq) {
      // std::cerr << "Reweighing " << (*it) << std::endl;
        Node *prev = sequences_[*it];
        Node *next = nullptr;
        while(prev) {
          next = nullptr;
          for(const auto& it2 : prev->outedges) {
            auto jt = std::find(it2->labels.begin(), it2->labels.end(), *it);
            if (jt != it2->labels.end()) {
              it2->weight += weight_adjust;
              next =it2->head;
              break;
            }
          }
          prev = next;
        }
      }
    }
  
    // Find the consensus and bundle
    TraverseHeaviestBundle();
    std::string dst{};
    for (const auto& it : consensus_) {
      dst += decoder_[it->code];
    }

    //std::cout << "Found concensus length is: " << dst.length() << "\n";

    // If consensus found
    if (dst.length() > 1){
      consensi->emplace_back(dst);  
    }
    
    // Remove the extra weights
    for(auto it = bundle.begin(); it != bundle.end(); ++it) {
      // std::cerr << "Reweighing " << (*it) << std::endl;
      Node *prev = sequences_[*it];
      Node *next = nullptr;
      while(prev) {
        next = nullptr;
        for(const auto& it2 : prev->outedges) {
          auto jt = std::find(it2->labels.begin(), it2->labels.end(), *it);
          if (jt != it2->labels.end()) {
            it2->weight -= weight_adjust;
            next =it2->head;
            break;
          }
        }
        prev = next;
      }
    }
  
    bundle.clear();

    // If consensus found
    if (dst.length() > 1){
      for(int it = 0; it < sequences_.size(); it++) {
        if (bundles[it] != 0){
          continue;
        }
        int current_score = 0;
        int good=0, all=0;
        for(auto n = sequences_[it]; n != nullptr; n = n->Successor(it)) {
          all++;
          auto ret = std::find(consensus_.begin(), consensus_.end(), n);
          if (ret != consensus_.end()) {
            good++;
          }
        }
        // std::cerr << good << " " << all << " " << (double)good/(double)all <<  std::endl;
        if ((double)good / (double)all >= bundling_threshold) {
          bundle.emplace_back(it);
          bundles[it] = bundle_id;
          num_bundled++;
        } else {
          if (it == best_seq){
            //std::cout << "Best seq not in its own bundle\n";
            bundles[it] = -1;
            num_bundled++;
          }
        }
      }
    } else {
      bundles[best_seq] = -1;
      num_bundled++;
    }

    // Set the weights to 0
    for(auto it = bundle.begin(); it != bundle.end(); ++it) {
      // std::cerr << "Reweighing " << (*it) << std::endl;
      Node *prev = sequences_[*it];
      Node *next = nullptr;
      while(prev) {
        next = nullptr;
        for(const auto& it2 : prev->outedges) {
          auto jt = std::find(it2->labels.begin(), it2->labels.end(), *it);
          if (jt != it2->labels.end()) {
          // ASSUMES ALL BASES HAD WEIGHT 1 INITIALLY
            it2->weight -= 2;  
            next =it2->head;
            break;
          }
        }
        prev = next;
      }
    }
    if (dst.length() > 1){
      bundle_id++;
      //std::cout << "Bundle id increased to " << bundle_id << "\n";
    }
  }
  delete [] seq_len;
}

  
void Graph::TraverseHeaviestBundle() {
  if (rank_to_node_.empty()) {
    return;
  }

  std::vector<Node*> predecessors(nodes_.size(), nullptr);
  std::vector<std::int64_t> scores(nodes_.size(), -1);
  Node* max = nullptr;

  for (const auto& it : rank_to_node_) {
    for (const auto& jt : it->inedges) {
      if ((scores[it->id] < jt->weight) ||
          (scores[it->id] == jt->weight && scores[predecessors[it->id]->id] <= scores[jt->tail->id])) {  // NOLINT
        scores[it->id] = jt->weight;
        predecessors[it->id] = jt->tail;
      }
    }
    if (predecessors[it->id]) {
      scores[it->id] += scores[predecessors[it->id]->id];
    }
    if (!max || scores[max->id] < scores[it->id]) {
      max = it;
    }
  }

  if (!max->outedges.empty()) {
    std::vector<std::uint32_t> node_id_to_rank(nodes_.size(), 0);
    for (std::uint32_t i = 0; i < rank_to_node_.size(); ++i) {
      node_id_to_rank[rank_to_node_[i]->id] = i;
    }
    while (!max->outedges.empty()) {
      max = BranchCompletion(node_id_to_rank[max->id], &scores, &predecessors);
    }
  }

  // traceback
  consensus_.clear();
  while (predecessors[max->id]) {
    consensus_.emplace_back(max);
    max = predecessors[max->id];
  }
  consensus_.emplace_back(max);
  std::reverse(consensus_.begin(), consensus_.end());
}

Graph::Node* Graph::BranchCompletion(
    std::uint32_t rank,
    std::vector<std::int64_t>* scores,
    std::vector<Node*>* predecessors) {
  auto start = rank_to_node_[rank];
  for (const auto& it : start->outedges) {
    for (const auto& jt : it->head->inedges) {
      if (jt->tail != start) {
        (*scores)[jt->tail->id] = -1;
      }
    }
  }

  Node* max = nullptr;
  for (std::uint32_t i = rank + 1; i < rank_to_node_.size(); ++i) {
    auto it = rank_to_node_[i];
    (*scores)[it->id] = -1;
    (*predecessors)[it->id] = nullptr;

    for (const auto& jt : it->inedges) {
      if ((*scores)[jt->tail->id] == -1) {
        continue;
      }
      if (((*scores)[it->id] < jt->weight) ||
          ((*scores)[it->id] == jt->weight && (*scores)[(*predecessors)[it->id]->id] <= (*scores)[jt->tail->id])) {  // NOLINT
        (*scores)[it->id] = jt->weight;
        (*predecessors)[it->id] = jt->tail;
      }
    }
    if ((*predecessors)[it->id]) {
      (*scores)[it->id] += (*scores)[(*predecessors)[it->id]->id];
    }
    if (!max || (*scores)[max->id] < (*scores)[it->id]) {
      max = it;
    }
  }

  return max;
}

std::vector<bool> Graph::ExtractSubgraph(const Node* begin, const Node* end) const {  // NOLINT
  std::vector<bool> dst(nodes_.size(), false);
  std::stack<const Node*> stack;
  stack.push(begin);

  while (!stack.empty()) {
    auto curr = stack.top();
    stack.pop();

    if (!dst[curr->id] && curr->id >= end->id) {
      for (const auto& it : curr->inedges) {
        stack.push(it->tail);
      }
      for (const auto& it : curr->aligned_nodes) {
        stack.push(it);
      }
      dst[curr->id] = true;
    }
  }

  return dst;
}

Graph Graph::Subgraph(
    std::uint32_t begin,
    std::uint32_t end,
    std::vector<const Node*>* subgraph_to_graph) const {
  if (!subgraph_to_graph) {
    throw std::invalid_argument(
        "[spoa::Graph::Subgraph] error: invalid ptr to subgraph_to_graph");
  }

  auto is_in_subgraph = ExtractSubgraph(nodes_[end].get(), nodes_[begin].get());

  // init subgraph
  Graph subgraph{};
  subgraph.num_codes_ = num_codes_;
  subgraph.coder_ = coder_;
  subgraph.decoder_ = decoder_;
  // subgraph.sequences_ = TODO(rvaser) maybe add sequences

  // create a map from subgraph nodes to graph nodes and vice versa
  subgraph_to_graph->clear();
  subgraph_to_graph->resize(nodes_.size(), nullptr);

  std::vector<Node*> graph_to_subgraph(nodes_.size(), nullptr);

  for (const auto& it : nodes_) {
    if (!is_in_subgraph[it->id]) {
      continue;
    }
    subgraph.AddNode(it->code);
    graph_to_subgraph[it->id] = subgraph.nodes_.back().get();
    (*subgraph_to_graph)[subgraph.nodes_.back()->id] = it.get();
  }

  // connect nodes
  for (const auto& it : nodes_) {
    if (!is_in_subgraph[it->id]) {
      continue;
    }
    auto jt = graph_to_subgraph[it->id];
    for (const auto& kt : it->inedges) {
      if (graph_to_subgraph[kt->tail->id]) {
        subgraph.AddEdge(graph_to_subgraph[kt->tail->id], jt, kt->weight);
      }
    }
    for (const auto& kt : it->aligned_nodes) {
      if (graph_to_subgraph[kt->id]) {
        jt->aligned_nodes.emplace_back(graph_to_subgraph[kt->id]);
      }
    }
  }

  subgraph.TopologicalSort();

  return subgraph;
}

void Graph::UpdateAlignment(
    const std::vector<const Node*>& subgraph_to_graph,
    Alignment* alignment) const {
  for (auto& it : *alignment) {
    if (it.first != -1) {
      it.first = subgraph_to_graph[it.first]->id;
    }
  }
}

void Graph::PrintDot(const std::string& path) const {
  if (path.empty()) {
    return;
  }
  std::ofstream os(path);

  std::vector<std::int32_t> consensus_rank(nodes_.size(), -1);
  std::int32_t rank = 0;
  for (const auto& it : consensus_) {
    consensus_rank[it->id] = rank++;
  }

  os << "digraph " << sequences_.size() << " {" << std::endl
     << "  graph [rankdir = LR]" << std::endl;
  for (const auto& it : nodes_) {
    os << "  " << it->id << "[label = \"" << it->id << " - "
       << static_cast<char>(decoder_[it->code]) << "\"";
    if (consensus_rank[it->id] != -1) {
      os << ", style = filled, fillcolor = goldenrod1";
    }
    os << "]" << std::endl;

    for (const auto& jt : it->outedges) {
      os << "  " << it->id << " -> " << jt->head->id
         << " [label = \"" << jt->weight << "\"";
      if (consensus_rank[it->id] + 1 == consensus_rank[jt->head->id]) {
        os << ", color = goldenrod1";
      }
      os << "]" << std::endl;
    }
    for (const auto& jt : it->aligned_nodes) {
      if (jt->id > it->id) {
        os << "  " << it->id << " -> " << jt->id
           << " [style = dotted, arrowhead = none]" << std::endl;
      }
    }
  }
  os << "}" << std::endl;

  os.close();
}

void Graph::Clear() {
  num_codes_ = 0;
  std::fill(coder_.begin(), coder_.end(), -1);
  std::fill(decoder_.begin(), decoder_.end(), -1);
  sequences_.clear();
  nodes_.clear();
  edges_.clear();
  rank_to_node_.clear();
  consensus_.clear();
}

}  // namespace spoa
