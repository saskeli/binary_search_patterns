#include <array>
#include <random>
#include <type_traits>
#include <vector>

#include "binary_search_patterns.hpp"
#include "counters/counters.hpp"

const constexpr size_t Q_COUNT = 100000;

template <class T, class searcher_t, size_t n, class counter_t>
void run_bench(counter_t& counter) {
  std::mt19937_64 gen;
  std::vector<T> arr;
  std::uniform_int_distribution<T> dist(0, std::numeric_limits<T>::max());
  for (size_t i = 0; i < n; ++i) {
    arr.push_back(dist(gen));
  }
  std::sort(arr.begin(), arr.end());
  searcher_t searcher(arr.data(), n);
  std::uniform_int_distribution<T> q_dist(arr[0],
                                          std::numeric_limits<T>::max() - 1);
  std::array<T, Q_COUNT> q_arr;
  for (size_t i = 0; i < q_arr.size(); ++i) {
    q_arr[i] = q_dist(gen);
  }
  uint64_t checksum = 0;
  counter.clear();
  for (auto q : q_arr) {
#ifdef DEPENDENCE_INSERTION
    q += checksum & 0b1;
#endif
    checksum += searcher.find(q);
  }
  counter.accumulate();
  counter.output_counters(0, Q_COUNT);
  std::cout << "Checsum: " << checksum << std::endl;
  std::cout << "Size (bytes): " << searcher.bytes() << std::endl;
}

#define TYPED_BENCH(typ, dtyp, n, block_size, short_circuit, prefetch)    \
  std::cout                                                               \
      << "\n====================== Profiling run ====================\n"; \
  std::cout << "Type: " << #typ << "\nData type: " << #dtyp               \
            << "\nSize: " << #n << "\nBlock_size: " << #block_size        \
            << "\nShort circuit: " << #short_circuit                      \
            << "\nPrefetch: " << #prefetch << std::endl;                  \
  run_bench<dtyp, typ<dtyp, block_size, short_circuit, prefetch>, n,      \
            decltype(counter)>(counter);

#define SS_BENCH(typ, dtyp, n, block_size, prefetch)    \
  TYPED_BENCH(typ, dtyp, n, block_size, true, prefetch) \
  TYPED_BENCH(typ, dtyp, n, block_size, false, prefetch)

#define SIZING_BENCH(typ, dtyp, block_size, prefetch) \
  SS_BENCH(typ, dtyp, 10000, block_size, prefetch)    \
  SS_BENCH(typ, dtyp, 100000, block_size, prefetch)   \
  SS_BENCH(typ, dtyp, 1000000, block_size, prefetch)  \
  SS_BENCH(typ, dtyp, 10000000, block_size, prefetch) \
  SS_BENCH(typ, dtyp, 100000000, block_size, prefetch)

#define DTYP_BENCH(typ, block_size, prefetch)       \
  SIZING_BENCH(typ, uint32_t, block_size, prefetch) \
  SIZING_BENCH(typ, uint64_t, block_size, prefetch)

#define PF_BENCH(typ, block_size)   \
  DTYP_BENCH(typ, block_size, true) \
  DTYP_BENCH(typ, block_size, false)

#define BLOCKING_BENCH(typ)    \
  PF_BENCH(typ, 4)             \
  PF_BENCH(typ, 8)             \
  PF_BENCH(typ, 16)            \
  PF_BENCH(typ, 32)            \
  PF_BENCH(typ, 64)            \
  PF_BENCH(typ, 128)           \
  DTYP_BENCH(typ, 256, false)  \
  DTYP_BENCH(typ, 512, false)  \
  DTYP_BENCH(typ, 1024, false) \
  DTYP_BENCH(typ, 2048, false) \
  DTYP_BENCH(typ, 4096, false) \
  DTYP_BENCH(typ, 8192, false)

int main() {
  count::Counters<false, 1, count::Counter::instructions,
                  count::Counter::branches, count::Counter::branch_miss,
                  count::Counter::L1D_miss, count::Counter::IPC>
      counter;
  DTYP_BENCH(control_binary_search, 1, false)
  DTYP_BENCH(heap_order_search, 1, false)
  BLOCKING_BENCH(b_plus_blocks)
  BLOCKING_BENCH(b_plus_heap_search)
  BLOCKING_BENCH(b_blocks)
  BLOCKING_BENCH(b_heap_search)
  return 0;
}