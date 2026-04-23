#include <benchmark/benchmark.h>

#include <array>
#include <random>
#include <type_traits>
#include <vector>

#include "binary_search_patterns.hpp"

template <class T, class searcher_t, size_t n>
void run_bench(benchmark::State& state) {
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
  std::array<T, 100000> q_arr;
  bool checksum_set = false;
  for (auto _ : state) {
    state.PauseTiming();
    for (size_t i = 0; i < q_arr.size(); ++i) {
      q_arr[i] = q_dist(gen);
    }
    uint64_t checksum = 0;
    state.ResumeTiming();
    for (auto q : q_arr) {
#ifdef DEPENDENCE_INSERTION
      q += checksum & 0b1;
#endif
      checksum += searcher.find(q);
    }
    state.PauseTiming();
    if (not checksum_set) {
      state.SetLabel(std::to_string(checksum));
      checksum_set = true;
    }
    state.ResumeTiming();
  }
}

#define TYPED_BENCH(typ, dtyp, n, block_size, short_circuit, prefetch)         \
  void BM_##typ##n##dtyp##block_size##short_circuit##prefetch(                 \
      benchmark::State& state) {                                               \
    run_bench<dtyp, typ<dtyp, block_size, short_circuit, prefetch>, n>(state); \
  }                                                                            \
  BENCHMARK(BM_##typ##n##dtyp##block_size##short_circuit##prefetch);

#define SS_BENCH(typ, dtyp, n, block_size, prefetch)    \
  TYPED_BENCH(typ, dtyp, n, block_size, true, prefetch) \
  TYPED_BENCH(typ, dtyp, n, block_size, false, prefetch)

#define DTYP_BENCH(typ, block_size, prefetch)              \
  SS_BENCH(typ, uint32_t, 10000, block_size, prefetch)     \
  SS_BENCH(typ, uint64_t, 10000, block_size, prefetch)     \
  SS_BENCH(typ, uint32_t, 100000, block_size, prefetch)    \
  SS_BENCH(typ, uint64_t, 100000, block_size, prefetch)    \
  SS_BENCH(typ, uint32_t, 1000000, block_size, prefetch)   \
  SS_BENCH(typ, uint64_t, 1000000, block_size, prefetch)   \
  SS_BENCH(typ, uint32_t, 10000000, block_size, prefetch)  \
  SS_BENCH(typ, uint64_t, 10000000, block_size, prefetch)  \
  SS_BENCH(typ, uint32_t, 100000000, block_size, prefetch) \
  SS_BENCH(typ, uint64_t, 100000000, block_size, prefetch)

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
  DTYP_BENCH(typ, 256, false)    \
  DTYP_BENCH(typ, 512, false)    \
  DTYP_BENCH(typ, 1024, false) \
  DTYP_BENCH(typ, 2048, false) \
  DTYP_BENCH(typ, 4096, false) \
  DTYP_BENCH(typ, 8192, false)

DTYP_BENCH(control_binary_search, 1, false)
DTYP_BENCH(heap_order_search, 1, false)
BLOCKING_BENCH(b_plus_blocks)
BLOCKING_BENCH(b_plus_heap_search)
BLOCKING_BENCH(b_blocks)
BLOCKING_BENCH(b_heap_search)

BENCHMARK_MAIN();