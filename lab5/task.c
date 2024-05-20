#include <stdio.h>
#include <stdlib.h>

#define TASK_COUNT              10000
#define TOTAL_SUM_WEIGHT        47000000
#define PROCESS_COUNT           16

int main() {
    FILE *file = fopen("task.txt", "w");

    int min_weight = 2 * TOTAL_SUM_WEIGHT / (TASK_COUNT * (PROCESS_COUNT + 1));
    printf("%d\n\n", min_weight);

    for (int i = 0; i < TASK_COUNT; ++i) {
        fprintf(file, "%d %d\n", i + 1, min_weight * (i % PROCESS_COUNT + 1));
    }

    int summary_weigh_1 = 0;
    int task_weight_process[PROCESS_COUNT] = {0};
    for (int i = 0; i < TASK_COUNT; ++i) {
        task_weight_process[i % PROCESS_COUNT] += min_weight * (i % PROCESS_COUNT + 1);
        summary_weigh_1 += min_weight * (i % PROCESS_COUNT + 1);
    }

    int summary_weigh_2 = 0;
    for (int i = 0; i < 16; ++i) {
        printf("%d\n", task_weight_process[i]);
        summary_weigh_2 += task_weight_process[i];
    }

    printf("%d %d", summary_weigh_1, summary_weigh_2);
}