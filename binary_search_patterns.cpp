#include <cstdint>
#include <iostream>
#include <vector>
#include <cstring>

class control_binary_search {
    private:
    uint64_t const* data_;
    const uint64_t n_;
   public:
    control_binary_search(uint64_t const* data, uint64_t n) : data_(data), n_(n) {}

    uint64_t find(uint64_t q) const {
        uint64_t a = 0;
        uint64_t b = n_ - 1;
        while (a < b) {
            uint64_t m = (a + b + 1) / 2;
            if (data_[m] < q) {
                a = m;
            } else {
                b = m - 1;
            }
        }
        return a;
    }
};

template <uint64_t block_size = 1024>
class b_star_blocks {
   private:
    static_assert(__builtin_popcountll(block_size) == 1);

    template <uint64_t block_size>
    class node {
       public:
        typedef std::pair<uint64_t, uint64_t> item;
        item children[block_size];

        node() {
            std::fill_n(children, block_size, {~uint64_t(0), 0});
        }

        node(const node& other) {
            std::memcpy(children, other.children, block_size * sizeof(sitem));
        }

        uint64_t build(uint64_t i_pos, uint64_t to_add, uint64_t n, uint64_t n_children, uint64_t c_level, node* nd_arr, uint64_t data) {
            for (uint64_t i = 0; i < n_children) {
                nd_arr[i + i_pos] = {};
                children[i] = {0, i_pos + i};
            }
            i_pos += n_children;
            for (uint64_t i = 0; i < n_children - 1) {
                
            }
        }
    };

    node<block_size>* nodes;
   public:
    b_star_blocks(uint64_t* data, uint64_t n) : nodes() {
        uint64_t nn = n;
        uint64_t node_count = 1;
        uint64_t levels = 0;
        while (nn > block_size) {
            node_count += nn / block_size;
            node_count += nn % block_size ? 1 : 0;
            levels++;
            nn /= block_size;
        }
        nodes = (node<block_size>*)malloc(node_count * sizeof(node<block_size>));
        uint64_t i = 0;
        nodes[0] = {};
        nodes[0].build(1, n, n, nn, levels, nodes, data);
    }

    uint64_t find(uint64_t q);
};

class heap_order_search {
    public:
    heap_order_search(uint64_t* data, uint64_t n) {

    }

    uint64_t find(uint64_t q);
};

int main(int argc, char const *argv[]) {
    if (argc < 3) {
        return 1;
    }
    uint64_t n, m, q;
    std::cin >> n >> m;
    uint64_t* data = (uint64_t*)malloc(n * sizeof(uint64_t));

    control_binary_search bins(data, n);
    b_star_blocks b_s(data, n);
    heap_order_search h_s(data, n);

    for (uint64_t i = 0; i < m; i++) {
        std::cin >> q;
        uint64_t r_bs = b_s.find(q);

        uint64_t r_hs = h_s.find(q);

        uint64_t r = bins.find(data, n, q);
    }

    free(data);
    return 0;
}
