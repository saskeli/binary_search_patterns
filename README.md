# binary_search_patterns

Testing a few binary search-based implementation on static sorted  data.

## binary search

Straight up classical binary search, with the implied inefficient memory access pattern.

Is generally faster than expected.

## B-star tree

A B-tree is built based on the data. Branch selection in internal nodes is done with a templated power-of-2 -based bianry search.

Is pretty ok but very space inefficient and loses comparative perofmance as the data structure size grows.

## Heap ordered binary search

Input data is "reordered" in such a way that the middle element is first, then the middle of the left and right sides, followed by the middles of the left and right portions of the left side and so on. For each "node" at index `i`, the left child is at `i * 2 + 1` and the right child is a t `i * 2 + 2`

Uses up to double the space of the raw data due to the data not necessarily conforming to the linearized tree topology.

Is typically the fastest of the approaches as long as the input data is large enough.
