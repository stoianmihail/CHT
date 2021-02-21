"Hist-Tree: Those Who Ignore It Are Doomed to Learn" - Andrew Crotty, Brown University
====

## Build

```
mkdir -p build
cd build
cmake -DCMAKE_BUILD_TYPE=Release ..
make
./example
```

## Examples

Using ``cht::Build`` to index the sorted data:

```c++
// Create random keys.
vector<uint64_t> keys(1e6);
generate(keys.begin(), keys.end(), rand);
keys.push_back(424242);
sort(keys.begin(), keys.end());

// Build CHT
const unsigned numBins = 64; // each node will have 64 separate bins
const unsigned maxError = 32; // the error of the index
cht::CHT cht(64, maxError);
for (const auto& key : keys) cht.AddKey(key);
cht.Build();

// Search using CHT
cht::SearchBound bound = cht.GetSearchBound(424242);
cout << "The search key is in the range: ["
			<< bound.begin << ", " << bound.end << ")" << endl;
auto start = begin(keys) + bound.begin, last = begin(keys) + bound.end;
auto pos = std::lower_bound(start, last, 424242) - begin(keys);
assert(keys[pos] == 424242);
cout << "The key is at position: " << pos << endl;
```
