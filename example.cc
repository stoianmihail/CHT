#include <iostream>
#include <algorithm>
#include <vector>

#include "include/cht/cht.h"

template <class KeyType>
void CompactHistTreeExample(bool useCache = false) {
  // Create random keys.
  std::vector<KeyType> keys(1e6);
  generate(keys.begin(), keys.end(), rand);
  keys.push_back(424242);
  std::sort(keys.begin(), keys.end());

  // Build CHT
  const unsigned numBins = 64; // each node will have 64 separate bins
  const unsigned maxError = 32; // the error of the index
  cht::CHT<KeyType> cht(numBins, maxError);
  for (const auto& key : keys) cht.AddKey(key);
	cht.Build();
	
  // Search using CHT
  cht::SearchBound bound = cht.GetSearchBound(424242);
  std::cout << "The search key is in the range: ["
       << bound.begin << ", " << bound.end << ")" << std::endl;
  auto start = std::begin(keys) + bound.begin, last = std::begin(keys) + bound.end;
	auto pos = std::lower_bound(start, last, 424242) - begin(keys);
	assert(keys[pos] == 424242);
  std::cout << "The key is at position: " << pos << std::endl;
}

int main(int argc, char** argv) {
  CompactHistTreeExample<uint32_t>();
  CompactHistTreeExample<uint64_t>();
  CompactHistTreeExample<uint32_t>(true);
	CompactHistTreeExample<uint64_t>(true);
	return 0;
}
