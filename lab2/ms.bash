gcc native.c -lm -Wall -Werror -Wextra -O3 -o native.out
time ./native.out

gcc openMPv$1.c -lm -fopenmp -Wall -Werror -Wextra -O3 -o siov$1.out
time ./siov$1.out