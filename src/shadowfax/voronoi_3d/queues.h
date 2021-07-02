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

#define QUEUE_TYPE int
#include "generic_lifo_queue.h"

#define QUEUE_TYPE int3
#include "generic_fifo_queue.h"

#endif  // SWIFTSIM_QUEUES_H
