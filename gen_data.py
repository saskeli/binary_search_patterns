#!python

import sys
import random

def main(n, m):
    r = random.Random()
    print(n)
    print(m)
    val = 0
    for _ in range(n):
        print(val)
        val += r.randint(1, 10)
    for _ in range(m):
        print(r.randint(0, val))


if __name__ == "__main__":
    if len(sys.argv) < 3:
        print("Number of values and queries required")
    else:
        main(int(sys.argv[1]), int(sys.argv[2]))