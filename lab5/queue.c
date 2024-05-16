#include <stdio.h>
#include <stdlib.h>

#include "queue.h"

Task_Queue *task_queue_create(int capacity) {
    Task_Queue *task_queue = malloc(sizeof(Task_Queue));

    if (task_queue == NULL) {
        return NULL;
    }

    Task *task_list = malloc (sizeof(Task) * capacity);
    if (task_list == NULL) {
        return NULL;
    }

    task_queue->task_list = task_list;
    task_queue->capacity = capacity;
    task_queue->size = 0;
    task_queue->pop_index = 0;

    return task_queue;
}

bool is_empty(const Task_Queue *task_queue) {
    return task_queue->size == 0;
}

bool is_full(const Task_Queue *task_queue) {
    return task_queue->size == task_queue->capacity;
}

int push(Task_Queue *task_queue, Task task) {
    if (task_queue == NULL) {
        return ERROR;
    }

    if (is_full(task_queue)) {
        return ERROR;
    }

    int push_index = (task_queue->pop_index + task_queue->size) % task_queue->capacity;
    task_queue->task_list[push_index] = task;
    task_queue->size++;

    return SUCCESS;
}

int pop(Task_Queue *task_queue, Task *task) {
    if (task_queue == NULL) {
        return ERROR;
    }

    if (is_empty(task_queue)) {
        return ERROR;
    }

    *task = task_queue->task_list[task_queue->pop_index];
    task_queue->pop_index = (task_queue->pop_index + 1) % task_queue->capacity;
    task_queue->size--;
    
    return SUCCESS;
}

void destroy(Task_Queue **task_queue) {
    if (*task_queue == NULL) {
        return;
    }

    if ((*task_queue)->task_list == NULL) {
        return;
    }

    free((*task_queue)->task_list);
    free(*task_queue);

    *task_queue = NULL;
}
