#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <stdbool.h>
#include <unistd.h>
#include <mpi.h>
// #include <mpich/mpi.h>
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
Task_Queue *task_queue;

pthread_mutex_t mutex;
pthread_cond_t worker_cond;
pthread_cond_t receiver_cond;

static inline void init_tasks() {
    const int TOTAL_SUM_WEIGHT = 50000000;
    int min_weight = 2 * TOTAL_SUM_WEIGHT / (TASK_COUNT * (process_count + 1));
    int task_id = 1;

    for (int i = 0; i < TASK_COUNT; ++i) {
        Task task = {
            .id = task_id,
            .process_id = process_id,
            .weight = min_weight * (i % process_count + 1)
        };

        if (i % process_count == process_id) {
            task_queue_push(task_queue, task);
            task_id++;
            process_sum_weight += task.weight;
        }
    }
}

static inline void execute_tasks() {
    while (true) {
        Task task;

        pthread_mutex_lock(&mutex);
        if (task_queue_is_empty(task_queue)) {
            pthread_mutex_unlock(&mutex);
            break;
        }
        task_queue_pop(task_queue, &task);
        pthread_mutex_unlock(&mutex);

        printf(FBLUE"Worker %d executing task %d of process %d and weight %d\n"FNORM,
               process_id,
               task.id,
               task.process_id,
               task.weight);
        usleep(task.weight);
    }
}

void *worker_start(void *args) {
    init_tasks();

    MPI_Barrier(MPI_COMM_WORLD);

    while (true) {
        execute_tasks();
        pthread_mutex_lock(&mutex);
        while (task_queue_is_empty(task_queue) && !termination) {
            pthread_cond_signal(&receiver_cond);
            pthread_cond_wait(&worker_cond, &mutex);
        }
        pthread_mutex_unlock(&mutex);
    }
}

void *receiver_start() {

}
 
void *sender_start() {

}

int main(int argc, char **argv) {
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

    start_time = MPI_Wtime();
    pthread_create(&worker_thread, NULL, worker_start, NULL);
    pthread_create(&receiver_thread, NULL, receiver_start, NULL);
    pthread_create(&sender_thread, NULL, sender_start, NULL);

    pthread_join(worker_thread, NULL);
    pthread_join(receiver_thread, NULL);
    pthread_join(sender_thread, NULL);
    end_time = MPI_Wtime();

    MPI_Barrier(MPI_COMM_WORLD);
    printf(FGREEN"Summary weight %d: %lf\n"FNORM, process_id, process_sum_weight * 1E-6);
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
