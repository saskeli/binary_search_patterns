#include <algorithm>
#include <array>
#include <cstddef>
#include <cstdint>
#include <iostream>
#include <random>

#include "binary_search_patterns.hpp"

template <class arr_t>
void report(arr_t& arr, uint32_t v, std::string typ, uint32_t expected,
            uint32_t actual, uint16_t block_size, bool short_circuit) {
  std::cout << "Problem with " << typ << ":\n";
  std::cout << "find(" << v << ") = " << actual << ", instead of expected "
            << expected << "\n";
  std::cout << "With block size " << block_size << " and short_circuit: " << short_circuit << "\n";
  std::cout << "const uint16_t size = " << arr.size() << ";\n";
  std::cout << "std::array<uint32_t, size> arr = {";
  for (size_t i = 0; i < arr.size(); ++i) {
    if (i % 8 == 0) {
      std::cout << "\n";
    }
    std::cout << arr[i] << (i == arr.size() - 1 ? "\n};" : ", ");
  }
  std::cout << "uint32_t q = " << v << ";\n"
            << "uint32_t res = " << expected << ";" << std::endl;
  std::cout << std::endl;
  exit(1);
}

template <class T, uint16_t block_size, T size, bool short_circuit>
void run_test(std::mt19937_64& gen, T limit) {
  std::vector<T> arr(size);
  std::uniform_int_distribution<T> uniform_dist(0, limit);
  for (size_t i = 0; i < size; ++i) {
    arr[i] = uniform_dist(gen);
  }
  std::sort(arr.begin(), arr.end());
  std::uniform_int_distribution<T> q_dist(arr[0], limit);
  control_binary_search<T, 1, short_circuit> b_s(arr.data(), size);
  branchless_binary_search<T, 1, short_circuit> bb_s(arr.data(), size);
  heap_order_search<T, 1, short_circuit> h_s(arr.data(), size);
  b_plus_blocks<T, block_size, short_circuit> b_plus(arr.data(), size);
  b_plus_heap_search<T, block_size, short_circuit> b_p_heap(arr.data(), size);
  b_blocks<T, block_size, short_circuit> b_b(arr.data(), size);
  b_heap_search<T, block_size, short_circuit> b_hs(arr.data(), size);
  for (size_t i = 0; i < 10000; ++i) {
    T v = q_dist(gen);
    T b_s_r = b_s.find(v);
    T bb_s_r = bb_s.find(v);
    T h_s_r = h_s.find(v);
    T b_plus_r = b_plus.find(v);
    T b_p_heap_r = b_p_heap.find(v);
    T b_b_r = b_b.find(v);
    T b_hs_r = b_hs.find(v);
    if (b_s_r != bb_s_r) {
      report(arr, v, "branchless binary search", b_s_r, bb_s_r, block_size, short_circuit);
    }
    if (b_s_r != h_s_r) {
      report(arr, v, "Heap-ordered binary search", b_s_r, h_s_r, block_size, short_circuit);
    }
    if (b_s_r != b_plus_r) {
      report(arr, v, "B+ blocks", b_s_r, b_plus_r, block_size, short_circuit);
    }
    if (b_s_r != b_p_heap_r) {
      report(arr, v, "heap ordered B+ blocks", b_s_r, b_p_heap_r, block_size, short_circuit);
    }
    if (b_s_r != b_b_r) {
      report(arr, v, "B blocks", b_s_r, b_b_r, block_size, short_circuit);
    }
    if (b_s_r != b_hs_r) {
      report(arr, v, "heap ordered B blocks", b_s_r, b_hs_r, block_size, short_circuit);
    }
  }
}

int main() {
  std::mt19937_64 gen;
  size_t epoc = 0;
  while (true) {
    run_test<uint32_t, 4, 100, false>(gen, 10000);
    run_test<uint32_t, 8, 100, true>(gen, 10000);
    run_test<uint32_t, 16, 100, false>(gen, 10000);
    run_test<uint32_t, 32, 100, true>(gen, 10000);
    run_test<uint32_t, 64, 100, false>(gen, 10000);

    run_test<uint32_t, 4, 1000, true>(gen, 10000);
    run_test<uint32_t, 8, 1000, false>(gen, 10000);
    run_test<uint32_t, 16, 1000, true>(gen, 10000);
    run_test<uint32_t, 32, 1000, false>(gen, 10000);
    run_test<uint32_t, 64, 1000, true>(gen, 10000);

    run_test<uint32_t, 4, 10000, false>(gen, 100000);
    run_test<uint32_t, 8, 10000, true>(gen, 100000);
    run_test<uint32_t, 16, 10000, false>(gen, 100000);
    run_test<uint32_t, 32, 10000, true>(gen, 100000);
    run_test<uint32_t, 64, 10000, false>(gen, 100000);

    run_test<uint64_t, 4, 1000, true>(gen, 100000);
    run_test<uint64_t, 8, 1000, false>(gen, 100000);
    run_test<uint64_t, 16, 1000, true>(gen, 100000);
    run_test<uint64_t, 32, 1000, false>(gen, 100000);
    run_test<uint64_t, 64, 1000, true>(gen, 100000);

    run_test<uint64_t, 4, 10000, false>(gen, 100000);
    run_test<uint64_t, 8, 10000, true>(gen, 100000);
    run_test<uint64_t, 16, 10000, false>(gen, 100000);
    run_test<uint64_t, 32, 10000, true>(gen, 100000);
    run_test<uint64_t, 64, 10000, false>(gen, 100000);

    run_test<uint64_t, 4, 100000, true>(gen, 1000000);
    run_test<uint64_t, 8, 100000, false>(gen, 1000000);
    run_test<uint64_t, 16, 100000, true>(gen, 1000000);
    run_test<uint64_t, 32, 100000, false>(gen, 1000000);
    run_test<uint64_t, 64, 100000, true>(gen, 1000000);
    std::cout << ++epoc << "\r" << std::flush;
  }
  return 0;
}
