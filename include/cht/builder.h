#pragma once

#include <cassert>
#include <cmath>
#include <limits>

#include "common.h"
#include "cht.h"

namespace cht {

// Allows building a `CompactHistTree` in a single pass over sorted data.
template <class KeyType>
class Builder {
 public:
	// The cache-oblivious structure makes sense when the tree becomes deep (`numBins` or `maxError` become small)
  Builder(KeyType min_key, KeyType max_key, size_t num_bins, size_t max_error, bool use_cache = false, bool single_pass = false)
      : min_key_(min_key),
        max_key_(max_key),
        num_bins_(num_bins),
        log_num_bins_(computeLog(num_bins_)),
        max_error_(max_error),
        use_cache_(use_cache),
        single_pass_(single_pass),
        curr_num_keys_(0),
        prev_key_(min_key) {
					assert((num_bins_ & (num_bins_ - 1)) == 0);
					// Compute the logarithm in base 2 of the range
					auto lg = computeLog(max_key_ - min_key_, true);
					
					// And also the initial shift for the first node of the tree
					assert(lg >= log_num_bins_);
					shift_ = lg - log_num_bins_;
				}

  // Adds a key. Assumes that keys are stored in a dense array.
  void AddKey(KeyType key) {
    assert(key >= min_key_ && key <= max_key_);
    // Keys need to be monotonically increasing.
    assert(key >= prev_key_);
    
		if (!single_pass_)
			keys_.push_back(key);
		else
			IncrementTable(key);
			
    ++curr_num_keys_;
    prev_key_ = key;
  }

  // Finalizes the construction and returns a read-only `RadixSpline`.
  CompactHistTree<KeyType> Finalize() {
    // Last key needs to be equal to `max_key_`.
    assert((!curr_num_keys_) || (prev_key_ == max_key_));

		if (!single_pass_)
			BuildOffline();
		
		// Flatten it
		if (!use_cache_)
			Flatten();
		else
			CacheObliviousFlatten();
		
    return CompactHistTree<KeyType>(
        min_key_, max_key_, curr_num_keys_,
				num_bins_, log_num_bins_, max_error_, shift_,
        std::move(table_));
  }

 private:
	static constexpr unsigned Infinity = std::numeric_limits<unsigned>::max();
	static constexpr unsigned Leaf = (1u << 31);
	static constexpr unsigned Mask = Leaf - 1;
	
	// Range covered by a node, i.e. [l, r[
	using Range = std::pair<unsigned, unsigned>;
	
	// (Node level, smallest key in node)
	using Info = std::pair<unsigned, KeyType>;

	static unsigned computeLog(uint32_t n, bool round = false) {
		assert(n);
		return 31 - __builtin_clz(n) + (round ? ((n & (n - 1)) != 0) : 0);
	}

	static unsigned computeLog(uint64_t n, bool round = false) {
		assert(n);
		return 63 - __builtin_clzl(n) + (round ? ((n & (n - 1)) != 0) : 0);
	}
	 
  void IncrementTable(KeyType key) {
    const auto Insert = [&]() -> void {
			// Traverse the tree from root
			for (unsigned level = 0, nodeIndex = 0; (shift_ >= level * log_num_bins_); ++level) {
				const auto [level, lower] = tree_[nodeIndex].first;
				
				// Compute the width and the bin for this node
				unsigned width = shift_ - level * log_num_bins_;
				auto bin = (key - min_key_ - tree_[nodeIndex].first.second) >> width;
				
				// Did we already visit this node?
				if (tree_[nodeIndex].second[bin].first != Infinity) {
					nodeIndex = tree_[nodeIndex].second[bin].first;
					continue;
				}
				
				// Can we continue with the next level?
				if (shift_ >= (level + 1) * log_num_bins_) {
					// Then set the partial sums, which will remain unchanged for this particular node
					tree_[nodeIndex].second[bin].first = curr_num_keys_;
					
					// Create the new node
					std::vector<Range> newNode;
					newNode.assign(num_bins_, {Infinity, Infinity});					
					
					// Compute the lowest key and attach the new node to the bin
					const auto newLower = lower + bin * (1ull << (shift_ - level * log_num_bins_));
					tree_.push_back({{level + 1, newLower}, newNode});
										
					// Point to the new node
					tree_[nodeIndex].second[bin].second = tree_.size() - 1;
					nodeIndex = tree_.size() - 1;
				}
			}
		};
		
		if (!curr_num_keys_)
			tree_.push_back({{0, 0}, std::vector<Range>(num_bins_, {Infinity, Infinity})});
		Insert();
  }

	void BuildOffline() {
		// Init the node, which covers the range `curr` := [a, b[
		auto initNode = [&](unsigned nodeIndex, Range curr) -> void {
			// Compute `width` of the current node (2^`width` represents the range covered by a single bin) 
			std::optional<unsigned> currBin = std::nullopt;
			unsigned width = shift_ - tree_[nodeIndex].first.first * log_num_bins_;
			
			// And compute the bins
			for (unsigned index = curr.first; index != curr.second; ++index) {
				// Extract the bin of the current key
				auto bin = (keys_[index] - min_key_ - tree_[nodeIndex].first.second) >> width;
				
				// Is the first bin or a new one?
				if ((!currBin.has_value()) || (bin != currBin.value())) {
					// Iterate the bins which has not been touched and set for them an empty range
					for (unsigned iter = currBin.has_value() ? (currBin.value() + 1) : 0; iter != bin; ++iter) {
						tree_[nodeIndex].second[iter] = {index, index};
					}
					
					// Init the current bin
					tree_[nodeIndex].second[bin] = {index, index};
					currBin = bin;
				}
				
				// And increase the range of the current bin
				tree_[nodeIndex].second[bin].second++;
			}
			assert(tree_[nodeIndex].second[currBin.value()].second == curr.second);
		};
		
		// Init the first node
		tree_.push_back({{0, 0}, std::vector<Range>(num_bins_, {curr_num_keys_, curr_num_keys_})});
		initNode(0, {0, curr_num_keys_});

		// Run the BFS
		std::queue<unsigned> nodes;
		nodes.push(0);
		while (!nodes.empty()) {
			// Extract from the queue
			auto node = nodes.front();
			nodes.pop();
			
			// Consider each bin and decide whether we should split it
			auto level = tree_[node].first.first;
			auto lower = tree_[node].first.second;
			for (unsigned index = 0; index != num_bins_; ++index) {
				// Should we split further?
				if (tree_[node].second[index].second - tree_[node].second[index].first >= max_error_) {
					// Corner-case: is #keys > range? Then create a leaf (this can only happen for datasets with duplicates)
					auto size = tree_[node].second[index].second - tree_[node].second[index].first;
					if (size >= (1ull << (shift_ - level * log_num_bins_))) {
						tree_[node].second[index].first |= Leaf;
						continue;
					}

					// Alloc the next node
					std::vector<Range> newNode;
					newNode.assign(num_bins_, {tree_[node].second[index].second, tree_[node].second[index].second});
					
					// And add it to the tree
					auto newLower = lower + index * (1ull << (shift_ - level * log_num_bins_));
					tree_.push_back({{level + 1, newLower}, newNode});
					
					// Init it
					initNode(tree_.size() - 1, tree_[node].second[index]);
					
					// Reset this node (no leaf, pointer to child)
					tree_[node].second[index] = {size, tree_.size() - 1};
					
					// And push it into the queue
					nodes.push(tree_.size() - 1);
				} else {
					// Leaf
					tree_[node].second[index].first |= Leaf;
				}
			}
		}
	}
	
	// Flatten the layout of the tree
	void Flatten() {
		table_.resize(tree_.size() * num_bins_);
		for (unsigned index = 0, limit = tree_.size(); index != limit; ++index) {
			for (unsigned bin = 0; bin != num_bins_; ++bin) {
				// Leaf node?
				if (tree_[index].second[bin].first & Leaf) {
					// Set the partial sum
					table_[(index << log_num_bins_) + bin] = tree_[index].second[bin].first;
				} else {
					// Set the pointer
					table_[(index << log_num_bins_) + bin] = (tree_[index].second[bin].second << log_num_bins_); 
				}
			}
		}
	}
	
	// Flatten the layout of the tree, such that the final layout is cache-oblivious
	void CacheObliviousFlatten() {
		// Build the precendence graph between nodes
		assert(!tree_.empty());
		auto maxLevel = tree_.back().first.first;
		std::vector<std::vector<unsigned>> graph(tree_.size());
		for (unsigned index = 0, limit = tree_.size(); index != limit; ++index) {
			graph[index].reserve(num_bins_);
			for (unsigned bin = 0; bin != num_bins_; ++bin) {
				// No leaf?
				if ((tree_[index].second[bin].first & Leaf) == 0) {
					graph[index].push_back(tree_[index].second[bin].second);
				}
			}
		}
		
		// And now set the count of nodes in subtree (bottom-up)
		auto access = [&](unsigned vertex) -> unsigned { return vertex * (maxLevel + 1); };
		std::vector<std::pair<unsigned, unsigned>> helper(tree_.size() * (maxLevel + 1), {Infinity, 0});
		for (unsigned index = 0, limit = tree_.size(); index != limit; ++index) {
			auto vertex = limit - index - 1;
			
			// Add the vertex itself
			helper[access(vertex) + tree_[vertex].first.first] = {vertex, 1};
			
			// And all subtrees, if any
			for (auto v : graph[vertex]) {
				for (unsigned lvl = 0; lvl <= maxLevel; ++lvl) {
					helper[access(vertex) + lvl].first = std::min(helper[access(vertex) + lvl].first, helper[access(v) + lvl].first);
					helper[access(vertex) + lvl].second += helper[access(v) + lvl].second;
				}
			}
		}
		
		// Build the `order`, the cache-oblivious permutation
		unsigned tempSize = 0;
		std::vector<unsigned> order(tree_.size());
		
		// Fill levels in [`lh`, `uh`[
		std::function<void(unsigned, unsigned, unsigned)> fill = [&](unsigned rowIndex, unsigned lh, unsigned uh) {
			// Stop?
			if (uh - lh == 1) {
				order[rowIndex] = tempSize++;
				return;
			}
			
			// Leaf?
			if (graph[rowIndex].empty()) {
				order[rowIndex] = tempSize++;
				return;
			}
			
			// Find the split level
			assert(helper[access(rowIndex) + lh].second);
			auto splitLevel = uh;
			while ((splitLevel >= lh + 1) && (!helper[access(rowIndex) + splitLevel - 1].second))
				--splitLevel;
			if (splitLevel == lh) {
				order[rowIndex] = tempSize++;
				return;
			}
			splitLevel = (lh + splitLevel) / 2; 

			// Recursion
			fill(rowIndex, lh, splitLevel);
			auto begin = helper[access(rowIndex) + splitLevel].first, end = begin + helper[access(rowIndex) + splitLevel].second;
			for (unsigned ptr = begin; ptr != end; ++ptr) {
				fill(ptr, splitLevel, uh); 
			}
		};
		
		// Start filling `order`, which is a permutation of the nodes, s.t. the tree becomes cache-oblivious
		fill(0, 0, maxLevel + 1);
		
		// Flatten with `order`
		table_.resize(tree_.size() * num_bins_);
		for (unsigned index = 0, limit = tree_.size(); index != limit; ++index) {
			for (unsigned bin = 0; bin != num_bins_; ++bin) {
				// Leaf node?
				if (tree_[index].second[bin].first & Leaf) {
					// Set the partial sum
					table_[(order[index] << log_num_bins_) + bin] = tree_[index].second[bin].first;
				} else {
					// Set the pointer
					table_[(order[index] << log_num_bins_) + bin] = (order[tree_[index].second[bin].second] << log_num_bins_); 
				}
			}
		}
		
		// And clean
		helper.clear();
		order.clear();
	}
  
  const KeyType min_key_;
  const KeyType max_key_;
  const size_t num_bins_;
	const size_t log_num_bins_;
  const size_t max_error_;
	const bool use_cache_;
	const bool single_pass_;
	
	size_t curr_num_keys_;
  KeyType prev_key_;
	size_t shift_;
  
	std::vector<KeyType> keys_;
	std::vector<unsigned> table_;
	std::vector<std::pair<Info, std::vector<Range>>> tree_;
};

}  // namespace cht
