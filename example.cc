#include <iostream>
#include <algorithm>
#include <vector>

#include "include/cht/builder.h"
#include "include/cht/cht.h"

template <class KeyType>
void CompactHistTreeExample(bool useCache = false) {
  // Create random keys.
#if 1
  std::vector<KeyType> keys(1e6);
  generate(keys.begin(), keys.end(), rand);
#else
	std::vector<KeyType> keys = {2, 3, 6, 10, 13, 17, 20, 21, 25, 30, 32, 40, 41, 51, 56, 60, 67, 68, 120, 140, 200, 400};
#endif
	static KeyType myKey = 424242;
  keys.push_back(myKey);
  std::sort(keys.begin(), keys.end());
  // Build `CompactHistTree`
	KeyType min = keys.front();
  KeyType max = keys.back();
  const unsigned numBins = 4; // each node will have 64 separate bins
  const unsigned maxError = 4; // the error of the index
  cht::Builder<KeyType> chtb(min, max, numBins, maxError, useCache, true);
  for (const auto& key : keys) chtb.AddKey(key);
	cht::CompactHistTree<KeyType> cht = chtb.Finalize();
	
  // Search using `CompactHistTree`
  cht::SearchBound bound = cht.GetSearchBound(myKey);
  std::cout << "The search key is in the range: ["
       << bound.begin << ", " << bound.end << ")" << std::endl;
  auto start = std::begin(keys) + bound.begin, last = std::begin(keys) + bound.end;
	auto pos = std::lower_bound(start, last, myKey) - begin(keys);
	assert(keys[pos] == myKey);
  std::cout << "The key is at position: " << pos << std::endl;
}

int main(int argc, char** argv) {
  CompactHistTreeExample<uint32_t>();
#if 0
  CompactHistTreeExample<uint64_t>();
  CompactHistTreeExample<uint32_t>(true);
	CompactHistTreeExample<uint64_t>(true);
#endif
	return 0;
}
