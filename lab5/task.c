#include <stdio.h>
#include <stdlib.h>

#define TASK_COUNT              10000
#define TOTAL_SUM_WEIGHT        47000000
#define PROCESS_COUNT           16

int main() {
    FILE *file = fopen("task.txt", "w");

    int min_weight = 2 * TOTAL_SUM_WEIGHT / (TASK_COUNT * (PROCESS_COUNT + 1));
    printf("%d", min_weight);

    for (int i = 0; i < TASK_COUNT; ++i) {
        fprintf(file, "%d %d\n", i + 1, min_weight * (i % PROCESS_COUNT + 1));
    }
}