#include <iostream>
#include <algorithm>
#include <vector>
#include <fstream>

#include "include/cht/builder.h"
#include "include/cht/cht.h"

template <class KeyType>
void CompactHistTreeExample(bool useCache = false) {
  // Create random keys.
#if 0
  std::vector<KeyType> keys(1e6);
  generate(keys.begin(), keys.end(), rand);
#else
#if 0
	std::vector<KeyType> keys = {2, 3, 6, 10, 13, 17, 20, 21, 25, 30, 32, 40, 41, 51, 56, 60, 67, 68, 120, 140, 200, 400, 500, 510, 520, 620, 640, 780, 1000, 1024, 1067, 1890, 1920, 2020, 2040, 2057, 3078, 4096};
#else
	std::ifstream in("keys.out");
	KeyType key;
	std::vector<KeyType> keys;
	while (in >> key) {
		keys.push_back(key);
	}
	
	std::ofstream out("debug.log");
	for (unsigned index = 677490; index != 677493; ++index) {
		out << "ba: " << keys[index] << std::endl;
	}
	//std::vector<KeyType> keys = {87102645, 87130124, 1314963317, 1400695062, 1485780449, 1762408883, 1986255891, 2252131116, 3746601961};
#endif
#endif
	static KeyType myKey = 424242;
  keys.push_back(myKey);
  std::sort(keys.begin(), keys.end());
  // Build `CompactHistTree`
	KeyType min = keys.front();
  KeyType max = keys.back();
  const unsigned numBins = 128; // each node will have 64 separate bins
  const unsigned maxError = 2; // the error of the index
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
