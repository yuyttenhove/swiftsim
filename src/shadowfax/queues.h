//
// Created by yuyttenh on 02/07/2021.
//

/**
 * @file queues.h
 *
 * @brief Generates code for a int LIFO queue and an int3 FIFO queue
 */

#ifndef SWIFTSIM_QUEUES_H
#define SWIFTSIM_QUEUES_H

#include "tuples.h"

#define QUEUE_SAFETY_CHECKS

#define QUEUE_TYPE int
#include "shadowfax/queues/generic_lifo_queue.h"
#undef QUEUE_TYPE

#define QUEUE_TYPE int2
#include "shadowfax/queues/generic_lifo_queue.h"
#undef QUEUE_TYPE

#define QUEUE_TYPE int3
#include "shadowfax/queues/generic_fifo_queue.h"

#endif  // SWIFTSIM_QUEUES_H
