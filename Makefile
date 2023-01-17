CL = $(shell getconf LEVEL1_DCACHE_LINESIZE)

binary_search_patterns: binary_search_patterns.cpp
	g++ -std=c++2a -march=native -DCACHE_LINE=$(CL) -Wall -Wextra -Wshadow -pedantic -Ofast -o binary_search_patterns binary_search_patterns.cpp

debug: binary_search_patterns.cpp
	g++ -std=c++2a -march=native -DCACHE_LINE=$(CL) -Wall -Wextra -Wshadow -pedantic -g -o binary_search_patterns binary_search_patterns.cpp