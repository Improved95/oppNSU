# gcc native.c -lm -Wall -Werror -Wextra -O3 -o native.out
# time ./native.out

gcc openMPv$1.c -lm -fopenmp -O3 -o siov$1.out
time OMP_NUM_THREADS=$2 ./siov$1.out