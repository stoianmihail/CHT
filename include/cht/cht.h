#pragma once

#include <algorithm>
#include <optional>
#include <cassert>
#include <vector>
#include <queue>
#include "common.h"

namespace cht {

class CHT {
 public:
  CHT() = default;

	// The constructor
  CHT(size_t numBins, size_t maxError) : num_bins_(numBins), max_error_(maxError) {}

  // Alternative for constructor
  void Init(size_t numBins, size_t maxError) {
    num_bins_ = numBins;
    max_error_ = maxError;
  }

  void AddKey(uint64_t key) {
		keys.push_back(key);
	}
  
  // Build the tree
  void Build() {
		// #bins should be a power of 2
		assert((num_bins_ & (num_bins_ - 1)) == 0);
		assert(std::is_sorted(keys.begin(), keys.end()));
		
		num_keys_ = keys.size();
		min_key_ = keys.front();
		max_key_ = keys.back();
		log_num_bins_ = computeLog(num_bins_);

		// Compute the logarithm in base 2 of the range
		auto diff = max_key_ - min_key_;
		auto lg = computeLog(diff);
		if (diff & (diff - 1))
			lg++;
		
		// And the initial shift for the first node of the tree
		assert(lg >= computeLog(num_bins_));
		shift_ = lg - computeLog(num_bins_);

		// Build the actual Compact Hist-Tree
		BuildCHT();
		
		// Flatten it
		Flatten();
		
		// And clear the storage
		keys.clear();
		cht.clear();
	}
	
  // Returns a search bound [`begin`, `end`) around the estimated position.
  SearchBound GetSearchBound(const uint64_t key) const {
    const size_t begin = Lookup(key);
		// `end` is exclusive.
    const size_t end = (begin + max_error_ + 1 > num_keys_) ? num_keys_ : (begin + max_error_ + 1);
    return SearchBound{begin, end};
	}

  // Returns the size in bytes.
  size_t GetSize() const {
    return sizeof(*this) + table.size() * sizeof(unsigned);
  }

 private:
	static constexpr unsigned LEAF = (1u << 31);
	static constexpr unsigned MASK = LEAF - 1;
	using range = std::pair<unsigned, unsigned>;
	using info = std::pair<unsigned, uint64_t>;

	static unsigned computeLog(uint32_t n) {
		assert(n);
		return 31 - __builtin_clz(n);
	}

	static unsigned computeLog(uint64_t n) {
		assert(n);
		return 63 - __builtin_clzl(n);
	}

	void BuildCHT() {
		// Init the node, which covers the range `curr` := [a, b[
		auto initNode = [&](unsigned nodeIndex, range curr) -> void {
			// Compute `width` of the current node (2^`width` represents the range covered by a single bin) 
			std::optional<unsigned> currBin = std::nullopt;
			unsigned width = shift_ - cht[nodeIndex].first.first * log_num_bins_;
			
			// And compute the bins
			for (unsigned index = curr.first; index != curr.second; ++index) {
				// Extract the bin of the current key
				auto bin = (keys[index] - min_key_ - cht[nodeIndex].first.second) >> width;
				
				// Is the first bin or a new one?
				if ((!currBin.has_value()) || (bin != currBin.value())) {
					// Iterate the bins which has not been touched and set for them an empty range
					for (unsigned iter = currBin.has_value() ? (currBin.value() + 1) : 0; iter != bin; ++iter) {
						cht[nodeIndex].second[iter] = std::make_pair(index, index);
					}
					
					// Init the current bin
					cht[nodeIndex].second[bin] = std::make_pair(index, index);
					currBin = bin;
				}
				
				// And increase the range of the current bin
				cht[nodeIndex].second[bin].second++;
			}
			assert(cht[nodeIndex].second[currBin.value()].second == curr.second);
		};
		
		// Init the first node
		cht.push_back(std::make_pair(std::make_pair(0, 0), std::vector<range>(num_bins_, std::make_pair(num_keys_, num_keys_))));
		initNode(0, std::make_pair(0, num_keys_));

		// Run the BFS
		std::queue<unsigned> nodes;
		nodes.push(0);
		while (!nodes.empty()) {
			// Extract from the queue
			auto node = nodes.front();
			nodes.pop();
			
			// Consider each bin and decide whether we should split it
			auto level = cht[node].first.first;
			auto lower = cht[node].first.second;
			for (unsigned index = 0; index != num_bins_; ++index) {
				// Should we split further?
				if (cht[node].second[index].second - cht[node].second[index].first >= max_error_) {
					// Corner-case: is #keys > range? Then create a leaf (this can only happen for datasets with duplicates)
					auto size = cht[node].second[index].second - cht[node].second[index].first;
					if (size >= (1ull << (shift_ - level * computeLog(num_bins_)))) {
						cht[node].second[index].first |= LEAF;
						continue;
					}

					// Alloc the next node
					std::vector<range> newNode;
					newNode.assign(num_bins_, std::make_pair(cht[node].second[index].second, cht[node].second[index].second));
					
					// And add it to the tree
					auto newLower = lower + index * (1ull << (shift_ - level * computeLog(num_bins_)));
					cht.push_back(std::make_pair(std::make_pair(level + 1, newLower), newNode));
					
					// Init it
					initNode(cht.size() - 1, cht[node].second[index]);
					
					// Reset this node (no leaf, pointer to child)
					cht[node].second[index] = std::make_pair(size, cht.size() - 1);
					
					// And push it into the queue
					nodes.push(cht.size() - 1);
				} else {
					// Leaf
					cht[node].second[index].first |= LEAF;
				}
			}
		}
	}

	// Flatten the layout of the tree
	void Flatten() {
		table.resize(cht.size() * num_bins_);
		for (unsigned index = 0, limit = cht.size(); index != limit; ++index) {
			for (unsigned bin = 0; bin != num_bins_; ++bin) {
				// Leaf node?
				if (cht[index].second[bin].first & LEAF) {
					// Set the partial sum
					table[(index << log_num_bins_) + bin] = cht[index].second[bin].first;
				} else {
					// Set the pointer
					table[(index << log_num_bins_) + bin] = (cht[index].second[bin].second << log_num_bins_); 
				}
			}
		}
	}
	 
	// Lookup `key` in tree
  size_t Lookup(uint64_t key) const {
		// Edge cases
		if (key <= min_key_) return 0;
		if (key >= max_key_) return num_keys_;
		key -= min_key_;
		
		auto width = shift_;
		auto next = 0;
		do {
			// Get the bin
			auto bin = key >> width;
			next = table[next + bin];

			// Is it a leaf?
			if (next & LEAF)
				return next & MASK;
			
			// Prepare for the next level
			key -= bin << width;
			width -= log_num_bins_;
		} while (true);
  }
	 
	size_t num_keys_;
  uint64_t min_key_;
  uint64_t max_key_;
	size_t max_error_ = 0;
	size_t num_bins_ = 0;
  size_t log_num_bins_;
  size_t shift_;
	std::vector<uint64_t> keys;
	std::vector<unsigned> table;
	std::vector<std::pair<info, std::vector<range>>> cht;	
};

} // namespace cht
