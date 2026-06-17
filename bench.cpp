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
  arr.push_back(0);
  for (size_t i = 1; i < n; ++i) {
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

template <class T, class searcher_t, size_t n>
void run_fp(benchmark::State& state) {
  std::mt19937_64 gen;
  std::vector<T> arr;
  std::uniform_real_distribution<T> dist(0.00003, 1.7e38);
  arr.push_back(0.0);
  for (size_t i = 1; i < n; ++i) {
    arr.push_back(dist(gen));
  }
  std::sort(arr.begin(), arr.end());
  searcher_t searcher(arr.data(), n);
  std::array<T, 100000> q_arr;
  bool checksum_set = false;
  for (auto _ : state) {
    state.PauseTiming();
    for (size_t i = 0; i < q_arr.size(); ++i) {
      q_arr[i] = dist(gen);
    }
    T checksum = 0.0;
    state.ResumeTiming();
    for (auto q : q_arr) {
#ifdef DEPENDENCE_INSERTION
      reinterpret_cast<char*>(&q)[0] ^= reinterpret_cast<char*>(&checksum)[0] & 0b1;
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

template <class T>
struct Decimal {
  T whole;
  T partial;

  static Decimal max_val() {
    return {std::numeric_limits<T>::max(), std::numeric_limits<T>::max()};
  }

  bool operator<(const Decimal& rhs) const {
    if (whole == rhs.whole) {
      return partial < rhs.partial;
    }
    return whole < rhs.whole;
  }

  bool operator>(const Decimal& rhs) const { return rhs < *this; }

  bool operator==(const Decimal& rhs) const {
    return whole == rhs.whole && partial == rhs.partial;
  }

  bool operator<=(const Decimal& rhs) const {
    if (whole == rhs.whole) {
      return partial <= rhs.partial;
    }
    return whole < rhs.whole;
  }

  bool operator>=(const Decimal& rhs) const {
    if (whole == rhs.whole) {
      return partial >= rhs.partial;
    }
    return whole > rhs.whole;
  }
};

template <class T, class searcher_t, size_t n>
void run_dec(benchmark::State& state) {
  std::mt19937_64 gen;
  std::vector<Decimal<T>> arr;
  std::uniform_int_distribution<T> dist(0, std::numeric_limits<T>::max() - 1);
  arr.push_back({0, 0});
  for (size_t i = 1; i < n; ++i) {
    arr.push_back({dist(gen), dist(gen)});
  }
  std::sort(arr.begin(), arr.end());
  searcher_t searcher(arr.data(), n);
  std::array<Decimal<T>, 100000> q_arr;
  bool checksum_set = false;
  for (auto _ : state) {
    state.PauseTiming();
    for (size_t i = 0; i < q_arr.size(); ++i) {
      q_arr[i] = {dist(gen), dist(gen)};
    }
    uint64_t checksum = 0;
    state.ResumeTiming();
    for (auto q : q_arr) {
#ifdef DEPENDENCE_INSERTION
      q.partial ^= checksum & 0b1;
#endif
      checksum += searcher.find(q).whole;
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

#define TYPED_FP_BENCH(typ, dtyp, n, block_size, short_circuit, prefetch)   \
  void FP_##typ##n##dtyp##block_size##short_circuit##prefetch(              \
      benchmark::State& state) {                                            \
    run_fp<dtyp, typ<dtyp, block_size, short_circuit, prefetch>, n>(state); \
  }                                                                         \
  BENCHMARK(FP_##typ##n##dtyp##block_size##short_circuit##prefetch);

#define TYPED_DEC_BENCH(typ, dtyp, n, block_size, short_circuit, prefetch)     \
  void DEC_##typ##n##dtyp##block_size##short_circuit##prefetch(               \
      benchmark::State& state) {                                              \
    run_dec<dtyp, typ<Decimal<dtyp>, block_size, short_circuit, prefetch>, n>( \
        state);                                                               \
  }                                                                           \
  BENCHMARK(DEC_##typ##n##dtyp##block_size##short_circuit##prefetch);

#define SS_BENCH(typ, dtyp, n, block_size, prefetch)    \
  TYPED_BENCH(typ, dtyp, n, block_size, true, prefetch) \
  TYPED_BENCH(typ, dtyp, n, block_size, false, prefetch)

#define SS_FP_BENCH(typ, dtyp, n, block_size, prefetch)    \
  TYPED_FP_BENCH(typ, dtyp, n, block_size, true, prefetch) \
  TYPED_FP_BENCH(typ, dtyp, n, block_size, false, prefetch)

#define SS_DEC_BENCH(typ, dtyp, n, block_size, prefetch)    \
  TYPED_DEC_BENCH(typ, dtyp, n, block_size, true, prefetch) \
  TYPED_DEC_BENCH(typ, dtyp, n, block_size, false, prefetch)

#define DTYP_BENCH(typ, block_size, prefetch)                  \
  SS_BENCH(typ, uint32_t, 10000, block_size, prefetch)         \
  SS_BENCH(typ, uint64_t, 10000, block_size, prefetch)         \
  SS_BENCH(typ, uint32_t, 100000, block_size, prefetch)        \
  SS_BENCH(typ, uint64_t, 100000, block_size, prefetch)        \
  SS_BENCH(typ, uint32_t, 1000000, block_size, prefetch)       \
  SS_BENCH(typ, uint64_t, 1000000, block_size, prefetch)       \
  SS_BENCH(typ, uint32_t, 10000000, block_size, prefetch)      \
  SS_BENCH(typ, uint64_t, 10000000, block_size, prefetch)      \
  SS_BENCH(typ, uint32_t, 100000000, block_size, prefetch)     \
  SS_BENCH(typ, uint64_t, 100000000, block_size, prefetch)     \
  SS_FP_BENCH(typ, float, 10000, block_size, prefetch)         \
  SS_FP_BENCH(typ, double, 10000, block_size, prefetch)        \
  SS_FP_BENCH(typ, float, 100000, block_size, prefetch)        \
  SS_FP_BENCH(typ, double, 100000, block_size, prefetch)       \
  SS_FP_BENCH(typ, float, 1000000, block_size, prefetch)       \
  SS_FP_BENCH(typ, double, 1000000, block_size, prefetch)      \
  SS_FP_BENCH(typ, float, 10000000, block_size, prefetch)      \
  SS_FP_BENCH(typ, double, 10000000, block_size, prefetch)     \
  SS_FP_BENCH(typ, float, 100000000, block_size, prefetch)     \
  SS_FP_BENCH(typ, double, 100000000, block_size, prefetch)    \
  SS_DEC_BENCH(typ, uint16_t, 10000, block_size, prefetch)     \
  SS_DEC_BENCH(typ, uint32_t, 10000, block_size, prefetch)     \
  SS_DEC_BENCH(typ, uint64_t, 10000, block_size, prefetch)     \
  SS_DEC_BENCH(typ, uint16_t, 100000, block_size, prefetch)    \
  SS_DEC_BENCH(typ, uint32_t, 100000, block_size, prefetch)    \
  SS_DEC_BENCH(typ, uint64_t, 100000, block_size, prefetch)    \
  SS_DEC_BENCH(typ, uint16_t, 1000000, block_size, prefetch)   \
  SS_DEC_BENCH(typ, uint32_t, 1000000, block_size, prefetch)   \
  SS_DEC_BENCH(typ, uint64_t, 1000000, block_size, prefetch)   \
  SS_DEC_BENCH(typ, uint16_t, 10000000, block_size, prefetch)  \
  SS_DEC_BENCH(typ, uint32_t, 10000000, block_size, prefetch)  \
  SS_DEC_BENCH(typ, uint64_t, 10000000, block_size, prefetch)  \
  SS_DEC_BENCH(typ, uint16_t, 100000000, block_size, prefetch) \
  SS_DEC_BENCH(typ, uint32_t, 100000000, block_size, prefetch) \
  SS_DEC_BENCH(typ, uint64_t, 100000000, block_size, prefetch)
  

#define PF_BENCH(typ, block_size)   \
  DTYP_BENCH(typ, block_size, true) \
  DTYP_BENCH(typ, block_size, false)

#define BLOCKING_BENCH(typ)   \
  PF_BENCH(typ, 4)            \
  PF_BENCH(typ, 8)            \
  PF_BENCH(typ, 16)           \
  PF_BENCH(typ, 32)           \
  PF_BENCH(typ, 64)           \
  PF_BENCH(typ, 128)          \
  DTYP_BENCH(typ, 256, false) \
  DTYP_BENCH(typ, 512, false) \
  DTYP_BENCH(typ, 1024, false)

DTYP_BENCH(control_binary_search, 1, false)
DTYP_BENCH(branchless_binary_search, 1, false)
DTYP_BENCH(heap_order_search, 1, false)
DTYP_BENCH(heap_order_search, 1, true)
BLOCKING_BENCH(b_plus_blocks_linear)
BLOCKING_BENCH(b_plus_blocks_logarithmic)
BLOCKING_BENCH(b_blocks_linear)
BLOCKING_BENCH(b_blocks_logarithmic)
BLOCKING_BENCH(b_plus_heap_linear)
BLOCKING_BENCH(b_plus_heap_logarithmic)
BLOCKING_BENCH(b_heap_linear)
BLOCKING_BENCH(b_heap_logarithmic)

BENCHMARK_MAIN();