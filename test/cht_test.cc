#include "include/cht/cht.h"

#include <random>
#include <unordered_set>

#include "gtest/gtest.h"
#include "include/cht/builder.h"

const size_t kNumKeys = 1000;
// Number of iterations (seeds) of random positive and negative test cases.
const size_t kNumIterations = 10;
const size_t kNumRadixBits = 18;
const size_t kMaxError = 32;

namespace {

// *** Helper methods ***

template <class KeyType>
std::vector<KeyType> CreateDenseKeys() {
  std::vector<KeyType> keys;
  keys.reserve(kNumKeys);
  for (size_t i = 0; i < kNumKeys; ++i) keys.push_back(i);
  return keys;
}

template <class KeyType>
std::vector<KeyType> CreateUniqueRandomKeys(size_t seed) {
  std::unordered_set<KeyType> keys;
  keys.reserve(kNumKeys);
  std::mt19937 g(seed);
  std::uniform_int_distribution<KeyType> d(std::numeric_limits<KeyType>::min(),
                                           std::numeric_limits<KeyType>::max());
  while (keys.size() < kNumKeys) keys.insert(d(g));
  std::vector<KeyType> sorted_keys(keys.begin(), keys.end());
  std::sort(sorted_keys.begin(), sorted_keys.end());
  return sorted_keys;
}

// Creates lognormal distributed keys, possibly with duplicates.
template <class KeyType>
std::vector<KeyType> CreateSkewedKeys(size_t seed) {
  std::vector<KeyType> keys;
  keys.reserve(kNumKeys);

  // Generate lognormal values.
  std::mt19937 g(seed);
  std::lognormal_distribution<double> d(/*mean*/ 0, /*stddev=*/2);
  std::vector<double> lognormal_values;
  lognormal_values.reserve(kNumKeys);
  for (size_t i = 0; i < kNumKeys; ++i) lognormal_values.push_back(d(g));
  const auto min_max =
      std::minmax_element(lognormal_values.begin(), lognormal_values.end());
  const double min = *min_max.first;
  const double max = *min_max.second;
  const double diff = max - min;

  // Scale values to the entire `KeyType` domain.
  const auto domain =
      std::numeric_limits<KeyType>::max() - std::numeric_limits<KeyType>::min();
  for (size_t i = 0; i < kNumKeys; ++i) {
    const double ratio = (lognormal_values[i] - min) / diff;
    keys.push_back(ratio * domain);
  }

  std::sort(keys.begin(), keys.end());
  return keys;
}

template <class KeyType>
cht::CompactHistTree<KeyType> CreateCompactHistTree(const std::vector<KeyType>& keys) {
  auto min = std::numeric_limits<KeyType>::min();
  auto max = std::numeric_limits<KeyType>::max();
  if (keys.size() > 0) {
    min = keys.front();
    max = keys.back();
  }
  cht::Builder<KeyType> chtb(min, max, kNumRadixBits, kMaxError);
  for (const auto& key : keys) chtb.AddKey(key);
  return chtb.Finalize();
}

template <class KeyType>
bool BoundContains(const std::vector<KeyType>& keys, cht::SearchBound bound,
                   KeyType key) {
  const auto it = std::lower_bound(keys.begin() + bound.begin,
                                   keys.begin() + bound.end, key);
  if (it == keys.end()) return false;
  return *it == key;
}

// *** Tests ***

template <class T>
struct CompactHistTreeTest : public testing::Test {
  using KeyType = T;
};

using AllKeyTypes = testing::Types<uint32_t, uint64_t>;
TYPED_TEST_SUITE(CompactHistTreeTest, AllKeyTypes);

TYPED_TEST(CompactHistTreeTest, AddAndLookupDenseKeys) {
  using KeyType = typename TestFixture::KeyType;
  const auto keys = CreateDenseKeys<KeyType>();
  const auto rs = CreateCompactHistTree(keys);
  for (const auto& key : keys)
    EXPECT_TRUE(BoundContains(keys, rs.GetSearchBound(key), key))
        << "key: " << key;
}

}  // namespace
