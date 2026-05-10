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

template <class T, class searcher_t, size_t n, class counter_t>
void run_fp_bench(counter_t& counter) {
  std::mt19937_64 gen;
  std::vector<T> arr;
  std::uniform_real_distribution<T> dist(0.0, std::numeric_limits<T>::max());
  arr.push_back(0.0);
  for (size_t i = 1; i < n; ++i) {
    arr.push_back(dist(gen));
  }
  std::sort(arr.begin(), arr.end());
  searcher_t searcher(arr.data(), n);
  std::array<T, Q_COUNT> q_arr;
  for (size_t i = 0; i < q_arr.size(); ++i) {
    q_arr[i] = dist(gen);
  }
  uint64_t checksum = 0;
  counter.clear();
  for (auto q : q_arr) {
#ifdef DEPENDENCE_INSERTION
    reinterpret_cast<char*>(&q)[0] ^= checksum & 0b1;
#endif
    checksum += searcher.find(q);
  }
  counter.accumulate();
  counter.output_counters(0, Q_COUNT);
  std::cout << "Checsum: " << checksum << std::endl;
  std::cout << "Size (bytes): " << searcher.bytes() << std::endl;
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

template <class T, class searcher_t, size_t n, class counter_t>
void run_dec_bench(counter_t& counter) {
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
  for (size_t i = 0; i < q_arr.size(); ++i) {
    q_arr[i] = {dist(gen), dist(gen)};
  }
  uint64_t checksum = 0;
  counter.clear();
  for (auto q : q_arr) {
#ifdef DEPENDENCE_INSERTION
    q.partial ^= checksum & 0b1;
#endif
    checksum += searcher.find(q).whole;
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

#define TYPED_FP_BENCH(typ, dtyp, n, block_size, short_circuit, prefetch) \
  std::cout                                                               \
      << "\n====================== Profiling run ====================\n"; \
  std::cout << "Type: " << #typ << "\nData type: " << #dtyp               \
            << "\nSize: " << #n << "\nBlock_size: " << #block_size        \
            << "\nShort circuit: " << #short_circuit                      \
            << "\nPrefetch: " << #prefetch << std::endl;                  \
  run_fp_bench<dtyp, typ<dtyp, block_size, short_circuit, prefetch>, n,   \
               decltype(counter)>(counter);

#define TYPED_DEC_BENCH(typ, dtyp, n, block_size, short_circuit, prefetch)     \
  std::cout                                                                    \
      << "\n====================== Profiling run ====================\n";      \
  std::cout << "Type: " << #typ << "\nData type: Decimal<" << #dtyp                    \
            << ">\nSize: " << #n << "\nBlock_size: " << #block_size             \
            << "\nShort circuit: " << #short_circuit                           \
            << "\nPrefetch: " << #prefetch << std::endl;                       \
  run_dec_bench<dtyp, typ<Decimal<dtyp>, block_size, short_circuit, prefetch>, \
                n, decltype(counter)>(counter);

#define SS_BENCH(typ, dtyp, n, block_size, prefetch)    \
  TYPED_BENCH(typ, dtyp, n, block_size, true, prefetch) \
  TYPED_BENCH(typ, dtyp, n, block_size, false, prefetch)

#define SS_FP_BENCH(typ, dtyp, n, block_size, prefetch)    \
  TYPED_FP_BENCH(typ, dtyp, n, block_size, true, prefetch) \
  TYPED_FP_BENCH(typ, dtyp, n, block_size, false, prefetch)

#define SS_DEC_BENCH(typ, dtyp, n, block_size, prefetch)    \
  TYPED_DEC_BENCH(typ, dtyp, n, block_size, true, prefetch) \
  TYPED_DEC_BENCH(typ, dtyp, n, block_size, false, prefetch)

#define SIZING_BENCH(typ, dtyp, block_size, prefetch) \
  SS_BENCH(typ, dtyp, 10000, block_size, prefetch)    \
  SS_BENCH(typ, dtyp, 100000, block_size, prefetch)   \
  SS_BENCH(typ, dtyp, 1000000, block_size, prefetch)  \
  SS_BENCH(typ, dtyp, 10000000, block_size, prefetch) \
  SS_BENCH(typ, dtyp, 100000000, block_size, prefetch)

#define SIZING_FP_BENCH(typ, dtyp, block_size, prefetch) \
  SS_FP_BENCH(typ, dtyp, 10000, block_size, prefetch)    \
  SS_FP_BENCH(typ, dtyp, 100000, block_size, prefetch)    \
  SS_FP_BENCH(typ, dtyp, 1000000, block_size, prefetch)   \
  SS_FP_BENCH(typ, dtyp, 10000000, block_size, prefetch)  \
  SS_FP_BENCH(typ, dtyp, 100000000, block_size, prefetch)

#define SIZING_DEC_BENCH(typ, dtyp, block_size, prefetch) \
  SS_DEC_BENCH(typ, dtyp, 10000, block_size, prefetch)    \
  SS_DEC_BENCH(typ, dtyp, 100000, block_size, prefetch)   \
  SS_DEC_BENCH(typ, dtyp, 1000000, block_size, prefetch)  \
  SS_DEC_BENCH(typ, dtyp, 10000000, block_size, prefetch) \
  SS_DEC_BENCH(typ, dtyp, 100000000, block_size, prefetch)

#define DTYP_BENCH(typ, block_size, prefetch)          \
  SIZING_BENCH(typ, uint32_t, block_size, prefetch)    \
  SIZING_BENCH(typ, uint64_t, block_size, prefetch)    \
  SIZING_FP_BENCH(typ, float, block_size, prefetch)    \
  SIZING_FP_BENCH(typ, double, block_size, prefetch)   \
  SIZING_DEC_BENCH(typ, uint16_t, block_size, prefetch) \
  SIZING_DEC_BENCH(typ, uint32_t, block_size, prefetch) \
  SIZING_DEC_BENCH(typ, uint64_t, block_size, prefetch)

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

int main() {
  count::Counters<false, 1, count::Counter::instructions,
                  count::Counter::branches, count::Counter::branch_miss,
                  count::Counter::L1D_miss, count::Counter::IPC>
      counter;
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

  return 0;
}