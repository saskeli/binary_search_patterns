#include <array>
#include <chrono>
#include <cstddef>
#include <cstring>
#include <iomanip>
#include <iostream>
#include <random>
#include <vector>

#include "search_microbench/searchers.hpp"

#ifndef CACHE_LINE
// Apparently the most common cache line size is 64.
#define CACHE_LINE 64
#endif

template <class T = uint64_t, uint16_t block_size = 1,
          bool short_circuit = false, bool prefetch = false>
class control_binary_search {
 private:
  T const* data_;
  const size_t n_;

 public:
  control_binary_search(T const* data, size_t n) : data_(data), n_(n) {}

  T find(T q) const {
    size_t a = 0;
    size_t b = n_ - 1;
    while (a < b) {
      size_t m = (a + b + 1) / 2;
      if constexpr (short_circuit) {
        if (data_[m] == q) [[unlikely]] {
          return q;
        }
      }
      if (data_[m] > q) {
        b = m - 1;
      } else {
        a = m;
      }
    }
    return data_[a];
  }

  size_t bytes() { return sizeof(control_binary_search) + n_ * sizeof(T); }
};

template <class T, uint16_t block_size = 1024, bool short_circuit = false,
          bool prefetch = false>
class b_plus_blocks {
 private:
  static_assert(__builtin_popcount(block_size) == 1);
  static_assert(block_size <= 8192);
  static_assert(block_size >= 2);

  class node {
   public:
    std::array<T, block_size> elements;
    std::array<node*, block_size> children;

    node() {
      std::fill_n(elements, block_size, std::numeric_limits<T>::max());
      std::fill_n(children, block_size, nullptr);
    }

    node(T const* arr, size_t elem_count = block_size) {
      std::copy_n(arr, elem_count, elements.data());
      if (elem_count < block_size) [[unlikely]] {
        std::fill_n(elements.data() + elem_count, block_size - elem_count,
                    std::numeric_limits<T>::max());
      }
      for (size_t i = 0; i < block_size; ++i) {
        children[i] = nullptr;
      }
    }

    node(node** arr, size_t elem_count = block_size) {
      for (size_t i = 0; i < elem_count; ++i) {
        elements[i] = arr[i]->elements[0];
        children[i] = arr[i];
      }
      if (elem_count < block_size) [[unlikely]] {
        std::fill_n(elements.data() + elem_count, block_size - elem_count,
                    std::numeric_limits<T>::max());
        std::fill_n(children.data() + elem_count, block_size - elem_count,
                    nullptr);
      }
    }

    node(const node& other) {
      elements = other.elements;
      children = other.children;
    }

    node& operator=(const node& other) {
      elements = other.elements;
      children = other.children;
      return *this;
    }

    ~node() {
      for (auto c : children) {
        if (c != nullptr) {
          delete (c);
        }
      }
    }

    size_t find(T q) const {
      if constexpr (prefetch) {
        const constexpr size_t node_bytes = sizeof(node);
        for (size_t i = 0; i < block_size; ++i) {
          uint8_t const* n_data =
              reinterpret_cast<uint8_t const*>(children.data() + 1);
          for (size_t cl_idx = 0; cl_idx < node_bytes; cl_idx += CACHE_LINE) {
            __builtin_prefetch(n_data + cl_idx);
          }
        }
      }
      return search(elements, q);
    }

    void print(std::string indent = "") {
      std::cout << indent;
      for (size_t i = 0; i < block_size; i++) {
        if (elements[i] == std::numeric_limits<T>::max()) {
          break;
        }
        if (i > 0 && i < block_size) {
          std::cout << ", ";
          if (i % 8 == 0) {
            std::cout << "\n" << indent;
          }
        }
        std::cout << elements[i];
      }
      std::cout << "\n" << std::endl;
      for (size_t i = 0; i < block_size; i++) {
        if (children[i] == nullptr) {
          break;
        }
        children[i]->print(indent + "   ");
      }
    }
  };

  node* root_;
  size_t levels_;
  size_t node_count_;

 public:
  b_plus_blocks(T const* data, size_t n) : levels_(1) {
    std::vector<size_t> c_per_level;
    size_t nn = n;
    node_count_ = 1;
    while (nn > block_size) {
      nn = (nn + block_size - 1) / block_size;
      node_count_ += nn;
      levels_++;
    }
    std::vector<node*> a_q;
    size_t i;
    node* nd_ptr;
    for (i = 0; i + block_size < n + 1; i += block_size) {
      nd_ptr = new node(data + i);
      a_q.push_back(nd_ptr);
    }
    if (n > i) {
      nd_ptr = new node(data + i, n - i);
      a_q.push_back(nd_ptr);
    }
    std::vector<node*> b_q;
    while (a_q.size() > 1) {
      for (i = 0; i + block_size < a_q.size() + 1; i += block_size) {
        nd_ptr = new node(a_q.data() + i);
        b_q.push_back(nd_ptr);
      }
      if (a_q.size() > i) {
        nd_ptr = new node(a_q.data() + i, a_q.size() - i);
        b_q.push_back(nd_ptr);
      }
      std::swap(a_q, b_q);
      b_q.clear();
    }
    root_ = a_q[0];
  }

  ~b_plus_blocks() {
    if (root_ != nullptr) {
      delete (root_);
    }
  }

  T find(T q) const {
    node* nd = root_;
    size_t idx;
    for (size_t i = 0; i < levels_ - 1; i++) {
      idx = nd->find(q);
      if constexpr (short_circuit) {
        if (nd->elements[idx] == q) [[unlikely]] {
          return q;
        }
      }
      nd = nd->children[idx];
    }
    idx = nd->find(q);
    return nd->elements[idx];
  }

  void print() {
    std::cout << "B-star search tree with " << node_count_ << " blocks"
              << std::endl;
    root_->print();
  }

  size_t bytes() { return sizeof(b_plus_blocks) + node_count_ * sizeof(node); }
};

template <class T, uint16_t block_size, bool short_circuit = false,
          bool prefetch = false>
class b_blocks {
 private:
  static_assert(__builtin_popcount(block_size) == 1);
  static_assert(block_size <= 8192);
  static_assert(block_size >= 2);

  class node {
   public:
    std::array<T, block_size> elements;
    std::array<node*, block_size> children;

    node() {
      std::fill_n(elements.data(), block_size, std::numeric_limits<T>::max());
      std::fill_n(children.data(), block_size, nullptr);
    }

    node(const node& other) {
      elements = other.elements;
      children = other.children;
    }

    node& operator=(const node& other) {
      elements = other.elements;
      children = other.children;
      return *this;
    }

    ~node() {
      for (auto c : children) {
        if (c != nullptr) {
          delete (c);
        }
      }
    }

    size_t find(T q) const {
      if constexpr (prefetch) {
        const constexpr size_t node_bytes = sizeof(node);
        for (size_t i = 0; i < block_size; ++i) {
          uint8_t const* n_data =
              reinterpret_cast<uint8_t const*>(children.data() + i);
          for (size_t c_idx = 0; c_idx < node_bytes; c_idx += CACHE_LINE) {
            __builtin_prefetch(n_data + c_idx);
          }
        }
      }
      return search(elements, q);
    }

    void print(std::string indent = "") {
      std::cout << indent;
      for (size_t i = 0; i < block_size; i++) {
        if (elements[i] == std::numeric_limits<T>::max()) {
          break;
        }
        if (i > 0 && i < block_size) {
          std::cout << ", ";
          if (i % 8 == 0) {
            std::cout << "\n" << indent;
          }
        }
        std::cout << elements[i];
      }
      std::cout << "\n" << std::endl;
      for (size_t i = 0; i < block_size; i++) {
        if (children[i] == nullptr) {
          break;
        }
        children[i]->print(indent + "   ");
      }
    }
  };

  node* root_;
  size_t levels_;
  size_t node_count_;

  template <class vec_t>
  void static build(T const* arr, size_t size, node* nd,
                    const vec_t& level_counts, size_t level,
                    size_t& node_count) {
    if (size <= block_size) {
      std::copy_n(arr, size, nd->elements.data());
      return;
    }
    for (size_t i = 0; i < block_size; ++i) {
      size_t t_idx = i * level_counts[level - 1] + i;
      if (t_idx >= size) {
        break;
      }
      nd->elements[i] = arr[t_idx];
      ++node_count;
      nd->children[i] = new node();
      size_t rest = size - t_idx - 1;
      size_t t_size =
          level_counts[level - 1] < rest ? level_counts[level - 1] : rest;
      build(arr + t_idx + 1, t_size, nd->children[i], level_counts, level - 1,
            node_count);
    }
  }

 public:
  b_blocks(T const* arr, size_t size) {
    std::vector<size_t> level_counts = {block_size};
    while (level_counts.back() < size) {
      level_counts.push_back(level_counts.back() * block_size + block_size);
    }
    levels_ = level_counts.size();
    root_ = new node();
    node_count_ = 1;
    build(arr, size, root_, level_counts, level_counts.size() - 1, node_count_);
  }

  ~b_blocks() {
    if (root_ != nullptr) {
      delete (root_);
    }
  }

  T find(T q) const {
    size_t idx = root_->find(q);
    T best = root_->elements[idx];
    node* nd = root_->children[idx];
    if constexpr (short_circuit) {
      if (best == q) [[unlikely]] {
        return q;
      }
    }
    for (size_t i = 1; i < levels_; i++) {
      if (nd == nullptr) [[unlikely]] {
        return best;
      }
      idx = nd->find(q);
      if (idx >= block_size) [[unlikely]] {
        return best;
      }
      T cand = nd->elements[idx];
      if (cand > q) [[unlikely]] {
        return best;
      }
      if constexpr (short_circuit) {
        if (nd->elements[idx] == q) [[unlikely]] {
          return q;
        }
      }
      if (cand <= q && cand > best) {
        best = cand;
      }
      nd = nd->children[idx];
    }
    return best;
  }

  void print() {
    std::cout << "B-tree with " << node_count_ << " blocks" << std::endl;
    root_->print();
  }

  size_t bytes() { return sizeof(b_blocks) + node_count_ * sizeof(node); }
};

template <class T, uint16_t block_size, bool short_circuit = false,
          bool prefetch = false>
class heap_order_search {
 public:
  std::vector<T> items_;

  heap_order_search(T const* data, size_t n)
      : items_((1 << (64 - __builtin_clzll(n))) - 1,
               std::numeric_limits<T>::max()) {
    size_t m = n / 2;
    items_[0] = data[m];
    build(0, m, 1, data);
    build(m + 1, n, 2, data);
  }

  T find(T q) const {
    size_t idx = 0;
    T pred = 0;
    while (idx < items_.size()) {
      if constexpr (short_circuit) {
        if (items_[idx] == q) [[unlikely]] {
          return q;
        }
      }
      T val = items_[idx];
      idx = val > q ? idx * 2 + 1 : (pred = std::max(pred, val), idx * 2 + 2);
    }
    return pred;
  }

  void print() {
    std::cout << "Heap order search tree with " << items_.size() << " nodes";
    for (size_t i = 0; i < items_.size(); i++) {
      if (i % 8 == 0) {
        std::cout << std::endl;
      }
      std::cout << items_[i] << "\t";
    }
    if ((items_.size() - 1) % 8 != 0) {
      std::cout << std::endl;
    }
  }

  size_t bytes() {
    return sizeof(heap_order_search) + items_.size() * sizeof(T);
  }

 private:
  void build(size_t left, size_t right, size_t target_index, T const* data) {
    if (left >= right) {
      return;
    }
    size_t m = (left + right) / 2;
    items_[target_index] = data[m];
    build(left, m, target_index * 2 + 1, data);
    build(m + 1, right, target_index * 2 + 2, data);
  }
};

template <class T, uint16_t block_size = 64, bool short_circuit = false,
          bool prefetch = false>
class b_plus_heap_search {
 private:
  static_assert(__builtin_popcountll(block_size) == 1);
  static_assert(block_size <= 8192);
  static_assert(block_size >= 2);

  class node {
   public:
    std::array<T, block_size> elements;

    node() {
      std::fill_n(elements.data(), block_size, std::numeric_limits<T>::max());
    }

    node(const node& other) { elements = other.elements; }

    node& operator=(const node& other) {
      elements = other.elements;
      return *this;
    }

    size_t find(uint64_t q) const { return search(elements, q); }

    void print() {
      for (uint64_t i = 0; i < block_size; i++) {
        std::cout << elements[i] << (i + 1 < block_size ? ", " : "");
      }
      std::cout << std::endl;
    }
  };

  node* nodes_;
  size_t levels_;
  size_t node_count_;

 public:
  b_plus_heap_search(T const* data, size_t n) : levels_(1) {
    size_t nn = (n + block_size - 1) / block_size;
    size_t leaves = nn;
    while (nn > block_size) {
      nn = (nn + block_size - 1) / block_size;
      levels_++;
    }
    size_t n_lev = 1;
    size_t t_nodes = n_lev;
    for (size_t i = 1; i < levels_; i++) {
      n_lev *= block_size;
      t_nodes += n_lev;
    }
    node_count_ = t_nodes + leaves;
    nodes_ = new node[node_count_]();
    std::memcpy(reinterpret_cast<T*>(nodes_ + t_nodes), data, n * sizeof(T));
    for (uint64_t i = t_nodes - 1; i < t_nodes; i--) {
      for (uint64_t k = 0; k < block_size; k++) {
        uint64_t c_idx = i * block_size + k + 1;
        if (c_idx < node_count_) {
          nodes_[i].elements[k] = nodes_[i * block_size + k + 1].elements[0];
        }
      }
    }
  }

  ~b_plus_heap_search() { delete[] (nodes_); }

  T find(T q) const {
    size_t n_idx = 0;
    for (size_t i = 0; i < levels_; i++) {
      if constexpr (prefetch) {
        const constexpr size_t tot_bytes = sizeof(node) * block_size;
        uint8_t const* byte_data =
            reinterpret_cast<uint8_t const*>(nodes_ + (n_idx * block_size + 1));
        for (size_t c_idx = 0; c_idx < tot_bytes; c_idx += CACHE_LINE) {
          __builtin_prefetch(byte_data + c_idx);
        }
      }
      uint16_t res = nodes_[n_idx].find(q);
      if constexpr (short_circuit) {
        if (nodes_[n_idx].elements[res] == q) {
          return q;
        }
      }
      n_idx = n_idx * block_size + 1 + nodes_[n_idx].find(q);
    }
    return nodes_[n_idx].elements[nodes_[n_idx].find(q)];
  }

  void print() {
    std::cout << "B-star heap search tree with " << node_count_ << " blocks"
              << std::endl;
    for (uint64_t i = 0; i < node_count_; i++) {
      std::cout << " node " << i << std::endl;
      nodes_[i].print();
    }
  }

  uint64_t bytes() {
    return sizeof(b_plus_heap_search) + node_count_ * sizeof(node);
  }
};

template <class T, uint16_t block_size = 64, bool short_circuit = false,
          bool prefetch = false>
class b_heap_search {
 private:
  static_assert(__builtin_popcountll(block_size) == 1);
  static_assert(block_size <= 8192);
  static_assert(block_size >= 2);

  class node {
   public:
    std::array<T, block_size> elements;

    node() {
      std::fill_n(elements.data(), block_size, std::numeric_limits<T>::max());
    }

    node(const node& other) { elements = other.elements; }

    node& operator=(const node& other) {
      elements = other.elements;
      return *this;
    }

    size_t find(uint64_t q) const { return search(elements, q); }

    void print() {
      for (uint64_t i = 0; i < block_size; i++) {
        std::cout << elements[i] << (i + 1 < block_size ? ", " : "");
      }
      std::cout << std::endl;
    }
  };

  std::vector<node> nodes_;
  size_t levels_;
  size_t node_count_;

  template <class vec_t, class lev_t>
  void static build(T const* arr, size_t size, vec_t& nodes, size_t idx,
                    const lev_t& level_counts, size_t level) {
    if (size <= block_size) {
      std::copy_n(arr, size, nodes[idx].elements.data());
      return;
    }
    for (size_t i = 0; i < block_size; ++i) {
      size_t t_idx = i * level_counts[level - 1] + i;
      if (t_idx >= size) {
        break;
      }
      nodes[idx].elements[i] = arr[t_idx];
      size_t child_idx = idx * block_size + 1 + i;
      if (nodes.size() <= child_idx) {
        nodes.resize(child_idx + 1);
      }
      size_t rest = size - t_idx - 1;
      size_t t_size =
          level_counts[level - 1] < rest ? level_counts[level - 1] : rest;
      build(arr + t_idx + 1, t_size, nodes, child_idx, level_counts, level - 1);
    }
  }

 public:
  b_heap_search(T const* arr, size_t size) : levels_(1) {
    std::vector<size_t> level_counts = {block_size};
    while (level_counts.back() < size) {
      level_counts.push_back(level_counts.back() * block_size + block_size);
    }
    levels_ = level_counts.size();
    nodes_ = {{}};
    build(arr, size, nodes_, 0, level_counts, level_counts.size() - 1);
    node_count_ = nodes_.size();
  }

  T find(T q) const {
    size_t idx = nodes_[0].find(q);
    T best = nodes_[0].elements[idx];
    size_t n_idx = 1 + idx;
    if constexpr (short_circuit) {
      if (best == q) {
        return q;
      }
    }
    for (size_t i = 1; i < levels_; i++) {
      if constexpr (prefetch) {
        uint8_t const* data_bytes = reinterpret_cast<uint8_t const*>(
            nodes_.data() + (n_idx * block_size + 1));
        const constexpr size_t tot_bytes = block_size * sizeof(node);
        for (size_t c_idx = 0; c_idx < tot_bytes; c_idx += CACHE_LINE) {
          __builtin_prefetch(data_bytes + c_idx);
        }
      }
      if (n_idx >= nodes_.size()) [[unlikely]] {
        return best;
      }
      idx = nodes_[n_idx].find(q);
      if (idx >= block_size) [[unlikely]] {
        return best;
      }
      T cand = nodes_[n_idx].elements[idx];
      if (cand > q) [[unlikely]] {
        return best;
      }
      if constexpr (short_circuit) {
        if (cand == q) [[unlikely]] {
          return q;
        }
      }
      if (cand <= q && cand > best) {
        best = cand;
      }
      n_idx = n_idx * block_size + 1 + idx;
    }
    return best;
  }

  void print() {
    std::cout << "B-heap search tree with " << node_count_ << " blocks"
              << std::endl;
    for (uint64_t i = 0; i < node_count_; i++) {
      std::cout << " node " << i << std::endl;
      nodes_[i].print();
    }
  }

  uint64_t bytes() {
    return sizeof(b_heap_search) + node_count_ * sizeof(node);
  }
};
