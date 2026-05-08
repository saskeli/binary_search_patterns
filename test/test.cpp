
#include <array>
#include <cstdint>
#include <string>

#include "../binary_search_patterns.hpp"
#include "../googletest/googletest/include/gtest/gtest.h"

template <class arr_t, uint16_t block_size, bool short_circuit>
void run_test(const arr_t& arr, typename arr_t::value_type q, typename arr_t::value_type res) {
  control_binary_search<typename arr_t::value_type, 1, short_circuit> b_s(arr.data(), arr.size());
  branchless_binary_search<typename arr_t::value_type, 1, short_circuit> bb_s(arr.data(), arr.size());
  heap_order_search<typename arr_t::value_type, 1, short_circuit> h_s(arr.data(), arr.size());
  b_plus_blocks<typename arr_t::value_type, block_size, short_circuit> b_star(arr.data(), arr.size());
  b_plus_heap_search<typename arr_t::value_type, block_size, short_circuit> b_s_heap(arr.data(), arr.size());
  b_blocks<typename arr_t::value_type, block_size, short_circuit> b_b(arr.data(), arr.size());
  b_heap_search<typename arr_t::value_type, block_size, short_circuit> b_hs(arr.data(), arr.size());

  ASSERT_EQ(b_s.find(q), res);
  ASSERT_EQ(bb_s.find(q), res);
  ASSERT_EQ(h_s.find(q), res);
  ASSERT_EQ(b_star.find(q), res);
  ASSERT_EQ(b_s_heap.find(q), res);
  ASSERT_EQ(b_b.find(q), res);
  ASSERT_EQ(b_hs.find(q), res);
}

TEST(Case, debug1) {
  const uint16_t size = 100;
  std::array<uint32_t, size> arr = {
      1927,  2207,  2271,  2903,  3091,  4553,  6983,  8070,  8641,  8886,
      9121,  11217, 13826, 14003, 15203, 16483, 19831, 21098, 22213, 22544,
      23927, 24553, 24916, 25048, 25131, 25635, 25726, 27296, 27419, 28411,
      28518, 29784, 29876, 32006, 32228, 34467, 36132, 36211, 36671, 40490,
      41937, 42295, 42381, 43414, 43809, 45564, 45992, 48096, 48683, 48733,
      49091, 49148, 49977, 50497, 51937, 52064, 52191, 52845, 53240, 54385,
      56103, 58388, 58866, 59698, 62041, 62380, 62890, 68782, 69623, 71067,
      71389, 71589, 72275, 73590, 73644, 74428, 75503, 75962, 76943, 77704,
      78682, 82830, 83636, 85321, 85707, 85990, 87039, 88879, 90711, 90924,
      91050, 92419, 93514, 94666, 95030, 96782, 96953, 98016, 98152, 98730};
  uint32_t q = 1928;
  uint32_t res = 1927;
  run_test<decltype(arr), 64, false>(arr, q, res);
  run_test<decltype(arr), 64, true>(arr, q, res);
}

TEST(Case, debug2) {
  const uint16_t size = 100;
  std::array<uint32_t, size> arr = {
      1927,  2207,  2271,  2903,  3091,  4553,  6983,  8070,  8641,  8886,
      9121,  11217, 13826, 14003, 15203, 16483, 19831, 21098, 22213, 22544,
      23927, 24553, 24916, 25048, 25131, 25635, 25726, 27296, 27419, 28411,
      28518, 29784, 29876, 32006, 32228, 34467, 36132, 36211, 36671, 40490,
      41937, 42295, 42381, 43414, 43809, 45564, 45992, 48096, 48683, 48733,
      49091, 49148, 49977, 50497, 51937, 52064, 52191, 52845, 53240, 54385,
      56103, 58388, 58866, 59698, 62041, 62380, 62890, 68782, 69623, 71067,
      71389, 71589, 72275, 73590, 73644, 74428, 75503, 75962, 76943, 77704,
      78682, 82830, 83636, 85321, 85707, 85990, 87039, 88879, 90711, 90924,
      91050, 92419, 93514, 94666, 95030, 96782, 96953, 98016, 98152, 98730};
  uint32_t q = 49091;
  uint32_t res = 49091;
  run_test<decltype(arr), 64, false>(arr, q, res);
  run_test<decltype(arr), 64, true>(arr, q, res);
}
