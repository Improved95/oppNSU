#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <stdbool.h>
#include <unistd.h>
#include <mpi.h>
#include <mpich/mpi.h>
#include <pthread.h>

#include "task_queue.h"
#include "color.h"

#define TASK_COUNT              2000
#define REQUEST_TAG             0
#define RESPONSE_TAG            1
#define EMPTY_QUEUE_RESPONSE    (-1)
#define TERMINATION_SIGNAL      (-2)

static int process_id, process_count;
static int process_sum_weight;
bool termination = false;
Task_Queue task_queue;

pthread_mutex_t mutex;
pthread_cond_t worker_cond;
pthread_cond_t receiver_cond;

int main() {
    int required = MPI_THREAD_MULTIPLE;
    int provided;
    double start_time;
    double end_time;
    pthread_t worker_thread;
    pthread_t receiver_thread;
    pthread_t sender_thread;

    MPI_Init_thread(&argc, &argv, required, &provided);
    if (provided != required) {
        return EXIT_FAILURE;
    }
    MPI_Comm_rank(MPI_COMM_WORLD, &process_id);
    MPI_Comm_size(MPI_COMM_WORLD, &process_count);

    task_queue = task_queue_create(TASK_COUNT);

    pthread_mutex_init(&mutex, NULL);
    pthread_cond_init(&worker_cond, NULL);
    pthread_cond_init(&receiver_cond, NULL);

    pthread_join(worker_thread, NULL);
    pthread_join(receiver_thread, NULL);
    pthread_join(sender_thread, NULL);
    end_time = MPI_Wtime();

    MPI_Barrier(MPI_COMM_WORLD);
    printf(FGREEN"Summary weight %d: %lf\n"FNORM, process_id, proc_sum_weight * 1E-6);
    MPI_Barrier(MPI_COMM_WORLD);
    if (process_id == 0) {
        printf(FGREEN"Time: %lf\n"FNORM, end_time - start_time);
    }

    task_queue_destroy(&task_queue);
    pthread_mutex_destroy(&mutex);
    pthread_cond_destroy(&worker_cond);
    pthread_cond_destroy(&receiver_cond);
    MPI_Finalize();

    return EXIT_SUCCESS;
}
