#include <iostream>
#include <chrono>
#include <random>

#include "include/cht/cht.h"

using namespace std::chrono;

template <class KeyType>
auto RandomNumberBetween = [](KeyType low, KeyType high) {
  auto randomFunc = [distribution_ = std::uniform_int_distribution<KeyType>(low, high), random_engine_ = std::mt19937{std::random_device{}()}]() mutable {
    return distribution_(random_engine_);
  };
  return randomFunc;
};

template <class KeyType>
static std::vector<KeyType> Generate(unsigned size, KeyType maxValue, bool withSort = false) {
  std::vector<KeyType> numbers;
  std::generate_n(std::back_inserter(numbers), size / 4, RandomNumberBetween<KeyType>(0, maxValue / 4));
  std::generate_n(std::back_inserter(numbers), size / 2, RandomNumberBetween<KeyType>(maxValue / 4, maxValue / 2));
  std::generate_n(std::back_inserter(numbers), size / 4, RandomNumberBetween<KeyType>(maxValue / 2, maxValue));
  if (withSort)
    sort(std::begin(numbers), std::end(numbers));
  return numbers;
}

template <class KeyType>
void CreateInput(std::vector<KeyType>& keys, std::vector<KeyType>& queries) {
	keys = Generate(5e6, std::numeric_limits<KeyType>::max(), true);
  queries = Generate(1e6, std::numeric_limits<KeyType>::max());
	std::cout << "Generated " << keys.size() << " keys and " << queries.size() << " lookup keys!" << std::endl; 
}

template <class KeyType>
void Benchmark() {
	std::vector<KeyType> keys, queries;
	CreateInput<KeyType>(keys, queries);
	
	for (unsigned i = 1; i <= 10; ++i) {
		for (unsigned j = 1; j <= 10; ++j) {
			auto numBins = 1u << i, maxError = 1u << j;
			
			// Build both CHTs
			cht::CHT<KeyType> cht(numBins, maxError, false);
			cht::CHT<KeyType> ccht(numBins, maxError, true);
			for (const auto& key : keys) cht.AddKey(key), ccht.AddKey(key); 
			cht.Build();
			ccht.Build();
			
			// Compare pure lookups
			auto measureTime = [&](std::string type) -> void {
				auto start = high_resolution_clock::now();
				if (type == "CHT") {
					for (auto query : queries) {
						cht.GetSearchBound(query);
					}
				} else if (type == "CCHT") {
					for (auto query : queries) {
						ccht.GetSearchBound(query);
					}
				}
				
				auto stop = high_resolution_clock::now();
				double answer = duration_cast<nanoseconds>(stop - start).count();
				std::cout << type << "<" << (std::is_same<KeyType, uint32_t>::value ? "uint32_t" : "uint64_t") << ">(numBins=" << numBins << ", maxError=" << maxError << "): " << answer << " ns" << std::endl;
			};
			
			for (unsigned index = 0; index != 5; ++index) {
				measureTime("CHT");
				measureTime("CCHT");
			}
		}
	}
}

int main(void) {
	Benchmark<uint32_t>();
	Benchmark<uint64_t>();
	return 0;
}
