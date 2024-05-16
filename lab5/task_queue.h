#ifndef TASK_QUEUE_H
#define TASK_QUEUE_H

#include <stdbool.h>

#define SUCCESS 0
#define ERROR   (-1)

typedef struct {
    int id;
    int process_id;
    int weight;
} Task;

typedef struct {
    Task *task_list;
    int capacity;
    int size;
    int pop_index;
} Task_Queue;

Task_Queue *task_queue_create(int capacity);
bool is_empty(const Task_Queue *task_queue);
bool is_full(const Task_Queue *task_queue);
int push(Task_Queue *task_queue, Task task);
int pop(Task_Queue *task_queue, Task *task);
void destroy(Task_Queue **task_queue);

#endif
