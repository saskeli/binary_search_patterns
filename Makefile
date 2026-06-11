.PHONY: clean run test

CL = $(shell getconf LEVEL1_DCACHE_LINESIZE)

MACHINE=$(shell lscpu | grep -o -P "(?<=Model name:).*" | sed -E 's/\s+//; s/\(\w+\)//g; s/\s/_/g')

CFLAGS=-march=native -std=c++2a -Wall -Wextra -Wshadow -pedantic
BENCH=-isystem benchmark/include -Lbenchmark/build/src -lbenchmark -lpthread
CMAKE_OPT=-DCMAKE_BUILD_TYPE=Release -DBENCHMARK_ENABLE_TESTING=OFF -DBENCHMARK_ENABLE_GTEST_TESTS=OFF -DBENCHMARK_DOWNLOAD_DEPENDENCIES=OFF

bench: bench.cpp binary_search_patterns.hpp search_microbench/searchers.hpp benchmark/build/src/libbenchmark.a
	g++ $(CFLAGS) -DCACHE_LINE=$(CL) -DNDEBUG -DDEPENDENCE_INSERTION -Ofast bench.cpp $(BENCH) -o bench

bench_avx: bench.cpp binary_search_patterns.hpp search_microbench/searchers.hpp benchmark/build/src/libbenchmark.a
	g++ $(CFLAGS) -DCACHE_LINE=$(CL) -DNDEBUG -Ofast bench.cpp $(BENCH) -o bench_avx

profile: profile.cpp binary_search_patterns.hpp search_microbench/searchers.hpp benchmark/build/src/libbenchmark.a
	g++ $(CFLAGS) -DCACHE_LINE=$(CL) -DNDEBUG -DDEPENDENCE_INSERTION -Ofast profile.cpp $(BENCH) -o profile

profile_avx: profile.cpp binary_search_patterns.hpp search_microbench/searchers.hpp benchmark/build/src/libbenchmark.a
	g++ $(CFLAGS) -DCACHE_LINE=$(CL) -DNDEBUG -Ofast profile.cpp $(BENCH) -o profile_avx

bf_test: bf_test.cpp binary_search_patterns.hpp search_microbench/searchers.hpp
	g++ $(CFLAGS) -DCACHE_LINE=$(CL) -DNDEBUG -Ofast bf_test.cpp -o bf_test

search_microbench/searchers.hpp:
	git submodule update --init

benchmark/include:
	git submodule update --init

benchmark/build/src/libbenchmark.a: | benchmark/include
	mkdir -p benchmark/build
	(cd benchmark; cmake $(CMAKE_OPT) -S . -B "build")
	(cd benchmark; cmake --build "build" --config Release)

googletest/googletest:
	git submodule update --init

googletest/build/lib/libgtest_main.a: | googletest/googletest
	(mkdir -p googletest/build && cd googletest/build && cmake .. && make)

test/test: googletest/build/lib/libgtest_main.a test/test.cpp binary_search_patterns.hpp search_microbench/searchers.hpp
	g++ -g $(CFLAGS) test/test.cpp -o test/test -lgtest_main -lgtest

test: test/test
	test/test $(ARG)

run: bench bench_avx profile profile_avx
	./bench | tee $(MACHINE).res
	./profile | tee $(MACHINE).prof

clean:
	rm -f bench
	rm -f bench_avx
	rm -f profile
	rm -f profile_avx