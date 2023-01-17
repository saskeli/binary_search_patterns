# Binary search patterns

Testing a few binary search-based implementation on static sorted  data.

### binary search

Straight up classical binary search, with the implied inefficient memory access pattern.

Is generally faster than expected.

### B-star tree

A B-tree is built based on the data. Branch selection in internal nodes is done with a templated power-of-2 -based binary search.

Is pretty ok but very space inefficient and loses comparative perofmance as the data structure size grows.

### Heap ordered binary search

Input data is "reordered" in such a way that the middle element is first, then the middle of the left and right sides, followed by the middles of the left and right portions of the left side and so on. For each "node" at index `i`, the left child is at `i * 2 + 1` and the right child is a t `i * 2 + 2`

Uses up to double the space of the raw data due to the data not necessarily conforming to the linearized tree topology.

Is typically the fastest of the approaches as long as the input data is large enough.

## Usage

Compile with ´make binary_search_patterns´.

Generate thest data with 

```bash
python gen_data.py <n> <m> > <file>
```

where `<n>` is the size of the data structure, ´<m>` the number of queries to generate and ´<file>´ and output file.

Run with

```bash
./binary_search_patterns <type> < <file> >> /dev/null
```

where ´<file>´ is a file generated by ´gen_data.py´ and ´<type>´ is

* 0 (or ommitted) to run queries interlieved on all the data structures at the same time,
* 1 to force cache invalidation between every query (**Caution: Very slow!)**,
* 2 to test each data structure separately.

Pipe data to ´/dev/null´ to avoid creating large tsv files with query data. If you want the per-query data, pipe somewhere else.
